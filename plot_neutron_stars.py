import gzip
import json
import argparse
import csv
from math import isnan
import time
import signal
import os

# Global flag for graceful shutdown
shutdown_requested = False


def signal_handler(sig, frame):
    global shutdown_requested
    shutdown_requested = True
    print("\n[SIGINT] Graceful shutdown requested...", flush=True)


def get_first(*keys, src=None):
    if not src:
        return None
    for k in keys:
        if k in src and src[k] is not None:
            return src[k]
    return None


def try_float(v):
    try:
        return float(v)
    except Exception:
        return None


def _clean_line_for_json(line):
    line = line.lstrip()
    if not line:
        return line
    if line == "[" or line == "]":
        return ""
    if line.startswith("["):
        line = line[1:].lstrip()
    while line.endswith(",") or line.endswith("]"):
        line = line[:-1].rstrip()
    return line


def extract_neutron_stars(gz_path, max_points=None, progress_interval=200000):
    global shutdown_requested
    points = []
    lines_read = 0
    systems_parsed = 0
    malformed_lines = 0
    stars_found = 0
    start = time.time()

    with gzip.open(gz_path, "rt", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            if shutdown_requested:
                elapsed = time.time() - start
                rate = lines_read / elapsed if elapsed > 0 else 0
                print(
                    f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] Shutdown requested — stopping read.",
                    flush=True,
                )
                print(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s rate={rate:.0f} L/s",
                    flush=True,
                )
                break

            lines_read += 1
            line = raw.strip()
            if not line:
                continue

            line = _clean_line_for_json(line)
            if not line:
                continue

            try:
                sysobj = json.loads(line)
            except Exception:
                malformed_lines += 1
            else:
                systems_parsed += 1
                sys_name = sysobj.get("name") or sysobj.get("systemName") or ""
                bodies = sysobj.get("bodies") or []
                for b in bodies:
                    if (
                        b.get("subType") != "Neutron Star"
                        and b.get("subtype") != "Neutron Star"
                    ):
                        continue
                    sr = try_float(
                        get_first("solarRadius", "radius", "stellarRadius", src=b)
                    )
                    rp = try_float(
                        get_first(
                            "rotationalPeriod", "rotationPeriod", "rotPeriod", src=b
                        )
                    )
                    age = try_float(
                        get_first("age", "age_Myr", "age_million_years", src=b)
                    )
                    mass = try_float(
                        get_first("solarMasses", "mass", "stellarMass", src=b)
                    )
                    temperature = try_float(
                        get_first("surfaceTemperature", "temperature", "temp", src=b)
                    )
                    if sr is None or rp is None or age is None:
                        continue
                    if isnan(sr) or isnan(rp) or isnan(age):
                        continue
                    if mass is not None and isnan(mass):
                        mass = None
                    if temperature is not None and isnan(temperature):
                        temperature = None

                    # Check for close companion that could contribute to accretion
                    has_close_companion = 0
                    parents = b.get("parents") or []
                    if parents and len(parents) > 0:
                        immediate_parent = parents[0]
                        # If immediate parent is Null (barycenter), it's in a binary system
                        if "Null" in immediate_parent:
                            barycenter_id = immediate_parent["Null"]
                            neutron_body_id = b.get("bodyId")
                            neutron_sma = try_float(
                                get_first("semiMajorAxis", "orbitalRadius", src=b)
                            )

                            # Find companions that share the same barycenter parent
                            for companion in bodies:
                                companion_parents = companion.get("parents") or []
                                if not companion_parents:
                                    continue
                                comp_immediate_parent = companion_parents[0]

                                # Check if this body shares the same barycenter
                                if "Null" in comp_immediate_parent:
                                    comp_barycenter_id = comp_immediate_parent["Null"]
                                    comp_body_id = companion.get("bodyId")

                                    if (
                                        comp_barycenter_id == barycenter_id
                                        and comp_body_id != neutron_body_id
                                    ):
                                        # Found a sibling in the binary system
                                        comp_type = companion.get("type") or ""
                                        comp_subtype = (
                                            companion.get("subType")
                                            or companion.get("subtype")
                                            or ""
                                        )

                                        # Only consider stellar companions (not planets)
                                        if comp_type == "Star":
                                            comp_sma = try_float(
                                                get_first(
                                                    "semiMajorAxis",
                                                    "orbitalRadius",
                                                    src=companion,
                                                )
                                            )
                                            comp_radius = try_float(
                                                get_first(
                                                    "solarRadius",
                                                    "radius",
                                                    "stellarRadius",
                                                    src=companion,
                                                )
                                            )

                                            # Calculate orbital separation
                                            if (
                                                neutron_sma is not None
                                                and comp_sma is not None
                                            ):
                                                # Total separation is sum of semi-major axes (in LS)
                                                separation_ls = neutron_sma + comp_sma
                                                # Convert to solar radii (1 LS ≈ 215.03 solar radii)
                                                separation_sr = separation_ls * 215.03

                                                # For accretion, typically need separation < ~1000 solar radii
                                                # and companion should be large enough (not a tiny star)
                                                # Conservative threshold: < 500 SR for close binaries
                                                if (
                                                    separation_sr < 500
                                                    and comp_radius is not None
                                                    and comp_radius > 0.1
                                                ):
                                                    has_close_companion = 1
                                                    break

                    points.append(
                        {
                            "system": sys_name,
                            "body": b.get("name") or b.get("bodyName") or "",
                            "solarRadius": sr,
                            "rotationalPeriod": rp,
                            "age": age,
                            "mass": mass,
                            "temperature": temperature,
                            "hasCloseCompanion": has_close_companion,
                        }
                    )
                    stars_found += 1
                    if max_points and len(points) >= max_points:
                        if progress_interval:
                            elapsed = time.time() - start
                            rate = lines_read / elapsed if elapsed > 0 else 0
                            print(
                                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s rate={rate:.0f} L/s",
                                flush=True,
                            )
                        return points

            if progress_interval and lines_read % progress_interval == 0:
                elapsed = time.time() - start
                rate = lines_read / elapsed if elapsed > 0 else 0
                print(
                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s rate={rate:.0f} L/s",
                    flush=True,
                )

    if progress_interval:
        elapsed = time.time() - start
        rate = lines_read / elapsed if elapsed > 0 else 0
        print(
            f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] DONE lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s avg_rate={rate:.0f} L/s",
            flush=True,
        )

    return points


def _range_for(vals, percentile_low=1, percentile_high=99):
    import numpy as np

    if len(vals) == 0:
        return (0.0, 1.0)
    arr = np.array(vals)
    mn = np.percentile(arr, percentile_low)
    mx = np.percentile(arr, percentile_high)
    if mn == mx:
        mn = arr.min()
        mx = arr.max()
        if mn == mx:
            mx = mn + (abs(mn) if mn != 0 else 1.0) * 1e-6
    return (mn, mx)


def write_csv_cache(points, csv_path):
    """Write neutron star data to CSV cache file"""
    start = time.time()
    print(f"Writing {len(points)} points to CSV: {csv_path}")

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        if not points:
            return csv_path

        # Get fieldnames from first point
        fieldnames = [
            "system",
            "body",
            "solarRadius",
            "rotationalPeriod",
            "age",
            "mass",
            "temperature",
            "hasCloseCompanion",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for point in points:
            writer.writerow(
                {
                    "system": point["system"],
                    "body": point["body"],
                    "solarRadius": point["solarRadius"],
                    "rotationalPeriod": point["rotationalPeriod"],
                    "age": point["age"],
                    "mass": point["mass"],
                    "temperature": point["temperature"],
                    "hasCloseCompanion": point["hasCloseCompanion"],
                }
            )

    elapsed = time.time() - start
    print(f"Wrote {len(points)} points to CSV in {elapsed:.1f}s")
    return csv_path


def read_csv_cache(csv_path):
    """Read neutron star data from CSV cache file"""
    start = time.time()
    print(f"Reading CSV cache from {csv_path}")

    points = []

    with open(csv_path, "r", encoding="utf-8") as f:
        csv_reader = csv.DictReader(f)
        for row in csv_reader:
            points.append(
                {
                    "system": row["system"],
                    "body": row["body"],
                    "solarRadius": float(row["solarRadius"]),
                    "rotationalPeriod": float(row["rotationalPeriod"]),
                    "age": float(row["age"]),
                    "mass": (
                        float(row["mass"])
                        if row["mass"] and row["mass"] != "None"
                        else None
                    ),
                    "temperature": (
                        float(row["temperature"])
                        if row["temperature"] and row["temperature"] != "None"
                        else None
                    ),
                    "hasCloseCompanion": int(row["hasCloseCompanion"]),
                }
            )

    elapsed = time.time() - start
    print(f"Read {len(points)} points from CSV in {elapsed:.1f}s")
    return points


def main():
    signal.signal(signal.SIGINT, signal_handler)

    p = argparse.ArgumentParser(description="Generate neutron star data CSV cache")
    p.add_argument(
        "--input",
        "-i",
        default="/home/meddler/spansh/galaxy.json.gz",
        help="path to galaxy.json.gz",
    )
    p.add_argument(
        "--output",
        "-o",
        default="neutron_stars_plots_cache.csv",
        help="output CSV file",
    )
    p.add_argument("--max", type=int, default=0, help="max points to read (0 = all)")
    p.add_argument(
        "--progress-interval",
        type=int,
        default=200000,
        help="lines between progress messages (0 to disable)",
    )

    args = p.parse_args()

    print("Extracting neutron stars...")
    pts = extract_neutron_stars(
        args.input,
        max_points=(args.max if args.max > 0 else None),
        progress_interval=args.progress_interval,
    )
    if not pts:
        print("No neutron star data found.")
        return

    point_count = len(pts)
    print(f"Extracted {point_count} neutron stars")

    # Write CSV cache
    print(f"Writing CSV to {args.output}...")
    write_csv_cache(pts, args.output)

    print(f"Done: {point_count} neutron stars cached to {args.output}")


if __name__ == "__main__":
    main()
