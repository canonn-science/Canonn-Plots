import gzip
import json
import argparse
import csv
from math import isnan
import time
import signal
import os
import sys
from collections import defaultdict

# Add EliteDangerousRegionMap to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "EliteDangerousRegionMap"))
from RegionMap import findRegion

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


def extract_neutron_stars(
    gz_path,
    systems_jsonl_path=None,
    max_points=None,
    progress_interval=200000,
    track_region_ages=False,
    region_ages_path="region_ages.json",
):
    global shutdown_requested
    points = []
    systems_written = 0
    lines_read = 0
    systems_parsed = 0
    malformed_lines = 0
    stars_found = 0
    start = time.time()

    # Initialize region age tracking
    region_stats = {}
    if track_region_ages:
        print(f"Region age tracking enabled, will save to: {region_ages_path}")
        # Initialize all regions 1-42
        for region_id in range(1, 43):
            region_stats[region_id] = {
                "id": region_id,
                "total_stars": 0,
                "stars_with_age": 0,
                "sum_age_myr": 0.0,
                "by_subtype": {},  # Will store {subType: {count, sum_age, stars_with_age}}
            }

    # Open temp file for writing systems incrementally
    systems_file = None
    systems_temp_path = None
    if systems_jsonl_path:
        systems_temp_path = systems_jsonl_path + ".tmp"
        systems_file = open(systems_temp_path, "w", encoding="utf-8")
        print(f"Writing systems to temporary file: {systems_temp_path}")

    try:
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

                    # Get system coordinates
                    coords = sysobj.get("coords") or {}
                    sys_x = coords.get("x")
                    sys_y = coords.get("y")
                    sys_z = coords.get("z")

                    system_has_neutron = False
                    bodies = sysobj.get("bodies") or []

                    # Track ALL stars for region age statistics
                    if (
                        track_region_ages
                        and sys_x is not None
                        and sys_y is not None
                        and sys_z is not None
                    ):
                        region_result = findRegion(sys_x, sys_y, sys_z)
                        if region_result:
                            region_id, _ = (
                                region_result  # Unpack tuple (id, region_data)
                            )
                            if 1 <= region_id <= 42:
                                for b in bodies:
                                    body_type = b.get("type") or ""
                                    if body_type == "Star":
                                        sub_type = (
                                            b.get("subType")
                                            or b.get("subtype")
                                            or "Unknown"
                                        )

                                        # Track overall stats
                                        region_stats[region_id]["total_stars"] += 1
                                        star_age = try_float(
                                            get_first(
                                                "age",
                                                "age_Myr",
                                                "age_million_years",
                                                src=b,
                                            )
                                        )
                                        if star_age is not None and not isnan(star_age):
                                            region_stats[region_id][
                                                "stars_with_age"
                                            ] += 1
                                            region_stats[region_id][
                                                "sum_age_myr"
                                            ] += star_age

                                        # Track per subType
                                        if (
                                            sub_type
                                            not in region_stats[region_id]["by_subtype"]
                                        ):
                                            region_stats[region_id]["by_subtype"][
                                                sub_type
                                            ] = {
                                                "total_count": 0,
                                                "stars_with_age": 0,
                                                "sum_age_myr": 0.0,
                                            }

                                        region_stats[region_id]["by_subtype"][sub_type][
                                            "total_count"
                                        ] += 1
                                        if star_age is not None and not isnan(star_age):
                                            region_stats[region_id]["by_subtype"][
                                                sub_type
                                            ]["stars_with_age"] += 1
                                            region_stats[region_id]["by_subtype"][
                                                sub_type
                                            ]["sum_age_myr"] += star_age

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
                            get_first(
                                "surfaceTemperature", "temperature", "temp", src=b
                            )
                        )
                        abs_magnitude = try_float(
                            get_first("absoluteMagnitude", "absMag", src=b)
                        )
                        if sr is None or rp is None or age is None:
                            continue
                        if isnan(sr) or isnan(rp) or isnan(age):
                            continue
                        if mass is not None and isnan(mass):
                            mass = None
                        if temperature is not None and isnan(temperature):
                            temperature = None
                        if abs_magnitude is not None and isnan(abs_magnitude):
                            abs_magnitude = None

                        # Determine spin direction based on environment
                        # Planets/rings cause spin-up, binary companions depend on mass ratio
                        spin = "none"
                        neutron_body_id = b.get("bodyId")

                        # Check for rings (indicate spin-up from accretion disk)
                        if b.get("rings") or b.get("hasRings"):
                            spin = "up"

                        # Check for planets orbiting the neutron star (spin-up from debris)
                        if spin == "none":
                            for body in bodies:
                                body_parents = body.get("parents") or []
                                body_type = body.get("type") or ""
                                # Check if this body is a planet orbiting our neutron star
                                if body_type == "Planet" and body_parents:
                                    first_parent = body_parents[0]
                                    # Check if neutron star is the direct parent
                                    if (
                                        "Star" in first_parent
                                        and first_parent["Star"] == neutron_body_id
                                    ):
                                        spin = "up"
                                        break

                        # Check for binary companion (mass-dependent spin evolution)
                        if spin == "none":
                            parents = b.get("parents") or []
                            if parents and len(parents) > 0:
                                immediate_parent = parents[0]
                                # If immediate parent is Null (barycenter), it's in a binary system
                                if "Null" in immediate_parent:
                                    barycenter_id = immediate_parent["Null"]
                                    neutron_sma = try_float(
                                        get_first(
                                            "semiMajorAxis", "orbitalRadius", src=b
                                        )
                                    )

                                    # Find stellar companions that share the same barycenter
                                    for companion in bodies:
                                        companion_parents = (
                                            companion.get("parents") or []
                                        )
                                        if not companion_parents:
                                            continue
                                        comp_immediate_parent = companion_parents[0]

                                        # Check if this body shares the same barycenter
                                        if "Null" in comp_immediate_parent:
                                            comp_barycenter_id = comp_immediate_parent[
                                                "Null"
                                            ]
                                            comp_body_id = companion.get("bodyId")

                                            if (
                                                comp_barycenter_id == barycenter_id
                                                and comp_body_id != neutron_body_id
                                            ):
                                                # Found a sibling in the binary system
                                                comp_type = companion.get("type") or ""

                                                # Only consider stellar companions
                                                if comp_type == "Star":
                                                    comp_sma = try_float(
                                                        get_first(
                                                            "semiMajorAxis",
                                                            "orbitalRadius",
                                                            src=companion,
                                                        )
                                                    )
                                                    comp_mass = try_float(
                                                        get_first(
                                                            "solarMasses",
                                                            "mass",
                                                            "stellarMass",
                                                            src=companion,
                                                        )
                                                    )

                                                    # Calculate orbital separation for close binary
                                                    if (
                                                        neutron_sma is not None
                                                        and comp_sma is not None
                                                    ):
                                                        # Total separation is sum of semi-major axes (in LS)
                                                        separation_ls = (
                                                            neutron_sma + comp_sma
                                                        )
                                                        # Convert to solar radii (1 LS ≈ 215.03 solar radii)
                                                        separation_sr = (
                                                            separation_ls * 215.03
                                                        )

                                                        # For accretion to occur, need close binary (< 500 SR)
                                                        if separation_sr < 500:
                                                            # Determine spin based on mass ratio
                                                            if (
                                                                mass is not None
                                                                and comp_mass
                                                                is not None
                                                            ):
                                                                # More massive companion → spin down
                                                                # (magnetic braking, complex torques)
                                                                # Less massive companion → spin up
                                                                # (standard recycling scenario)
                                                                if comp_mass > mass:
                                                                    spin = "down"
                                                                else:
                                                                    spin = "up"
                                                                break
                                                            elif (
                                                                comp_mass is not None
                                                                and comp_mass > 0.5
                                                            ):
                                                                # If we don't know neutron star mass,
                                                                # assume typical NS mass ~1.4 M☉
                                                                if comp_mass > 1.4:
                                                                    spin = "down"
                                                                else:
                                                                    spin = "up"
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
                                "absoluteMagnitude": abs_magnitude,
                                "x": sys_x,
                                "y": sys_y,
                                "z": sys_z,
                                "spin": spin,
                            }
                        )
                        stars_found += 1
                        system_has_neutron = True
                        if max_points and len(points) >= max_points:
                            if progress_interval:
                                elapsed = time.time() - start
                                rate = lines_read / elapsed if elapsed > 0 else 0
                                print(
                                    f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} malformed={malformed_lines} elapsed={elapsed:.0f}s rate={rate:.0f} L/s",
                                    flush=True,
                                )
                            return points

                # Write system to jsonl if it has neutron stars
                if (
                    system_has_neutron
                    and sys_x is not None
                    and sys_y is not None
                    and sys_z is not None
                ):
                    if systems_file:
                        # Write system to file immediately
                        systems_file.write(
                            json.dumps(
                                {
                                    "name": sys_name,
                                    "coords": {"x": sys_x, "y": sys_y, "z": sys_z},
                                    "bodies": bodies,
                                }
                            )
                            + "\n"
                        )
                        systems_written += 1

                if progress_interval and lines_read % progress_interval == 0:
                    elapsed = time.time() - start
                    rate = lines_read / elapsed if elapsed > 0 else 0
                    print(
                        f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} systems_saved={systems_written} malformed={malformed_lines} elapsed={elapsed:.0f}s rate={rate:.0f} L/s",
                        flush=True,
                    )

        if progress_interval:
            elapsed = time.time() - start
            rate = lines_read / elapsed if elapsed > 0 else 0
            print(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] DONE lines={lines_read} systems={systems_parsed} neutron_stars={stars_found} systems_saved={systems_written} malformed={malformed_lines} elapsed={elapsed:.0f}s avg_rate={rate:.0f} L/s",
                flush=True,
            )

    finally:
        # Close and rename temp file
        if systems_file:
            systems_file.close()
            if systems_temp_path and os.path.exists(systems_temp_path):
                # Rename temp file to final name
                if os.path.exists(systems_jsonl_path):
                    os.remove(systems_jsonl_path)
                os.rename(systems_temp_path, systems_jsonl_path)
                print(f"Saved {systems_written} systems to {systems_jsonl_path}")

        # Save region age statistics
        if track_region_ages:
            print(f"Calculating region averages and saving to {region_ages_path}...")

            # Load region names from RegionMap
            region_map_path = os.path.join(
                os.path.dirname(__file__),
                "EliteDangerousRegionMap",
                "RegionMapData.json",
            )
            region_names = {}
            try:
                with open(region_map_path, "r") as f:
                    region_data = json.load(f)
                    regions = region_data.get("regions", [])
                    for i, region in enumerate(regions):
                        if region and "name" in region:
                            region_names[i] = region["name"]
            except Exception as e:
                print(f"Warning: Could not load region names: {e}")

            # Collect all unique star subTypes across all regions
            all_subtypes = set()
            for region_id in range(1, 43):
                all_subtypes.update(region_stats[region_id]["by_subtype"].keys())

            # Build output structure
            output = {
                "metadata": {
                    "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
                    "total_systems_processed": systems_parsed,
                    "regions_included": 42,
                    "star_subtypes": sorted(list(all_subtypes)),
                },
                "regions": {},
            }

            total_stars_with_age = 0
            for region_id in range(1, 43):
                stats = region_stats[region_id]
                avg_age = 0.0
                if stats["stars_with_age"] > 0:
                    avg_age = stats["sum_age_myr"] / stats["stars_with_age"]
                    total_stars_with_age += stats["stars_with_age"]

                # Build subtype data
                subtypes_data = {}
                for subtype, subtype_stats in stats["by_subtype"].items():
                    avg_subtype_age = 0.0
                    if subtype_stats["stars_with_age"] > 0:
                        avg_subtype_age = (
                            subtype_stats["sum_age_myr"]
                            / subtype_stats["stars_with_age"]
                        )

                    subtypes_data[subtype] = {
                        "total_count": subtype_stats["total_count"],
                        "stars_with_age": subtype_stats["stars_with_age"],
                        "sum_age_myr": round(subtype_stats["sum_age_myr"], 2),
                        "average_age_myr": round(avg_subtype_age, 2),
                    }

                output["regions"][str(region_id)] = {
                    "id": region_id,
                    "name": region_names.get(region_id, f"Region {region_id}"),
                    "total_stars": stats["total_stars"],
                    "stars_with_age": stats["stars_with_age"],
                    "sum_age_myr": stats["sum_age_myr"],
                    "average_age_myr": round(avg_age, 2),
                    "by_subtype": subtypes_data,
                }

            output["metadata"]["total_stars_with_age"] = total_stars_with_age

            # Write to file
            with open(region_ages_path, "w") as f:
                json.dump(output, f, indent=2)

            print(f"Saved region age statistics to {region_ages_path}")
            print(f"Total stars with age data: {total_stars_with_age:,}")
            print(f"Unique star subTypes found: {len(all_subtypes)}")
            print(f"SubTypes: {', '.join(sorted(all_subtypes))}")

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


def write_systems_jsonl(systems, jsonl_path):
    """Write systems with neutron stars to JSONL file"""
    start = time.time()
    print(f"Writing {len(systems)} systems to JSONL: {jsonl_path}")

    with open(jsonl_path, "w", encoding="utf-8") as f:
        for system in systems:
            f.write(json.dumps(system) + "\n")

    elapsed = time.time() - start
    print(f"Wrote {len(systems)} systems to JSONL in {elapsed:.1f}s")
    return jsonl_path


def load_systems_jsonl(jsonl_path):
    """Load systems from JSONL file"""
    start = time.time()
    print(f"Loading systems from JSONL: {jsonl_path}")

    systems = []
    with open(jsonl_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                systems.append(json.loads(line))

    elapsed = time.time() - start
    print(f"Loaded {len(systems)} systems from JSONL in {elapsed:.1f}s")
    return systems


def extract_neutron_stars_from_systems(systems, max_points=None):
    """Extract neutron stars from pre-loaded system data"""
    print(f"Extracting neutron stars from {len(systems)} systems...")
    points = []
    stars_found = 0

    for sysobj in systems:
        sys_name = sysobj.get("name") or sysobj.get("systemName") or ""
        coords = sysobj.get("coords") or {}
        sys_x = coords.get("x")
        sys_y = coords.get("y")
        sys_z = coords.get("z")

        bodies = sysobj.get("bodies") or []
        for b in bodies:
            if (
                b.get("subType") != "Neutron Star"
                and b.get("subtype") != "Neutron Star"
            ):
                continue

            sr = try_float(get_first("solarRadius", "radius", "stellarRadius", src=b))
            rp = try_float(
                get_first("rotationalPeriod", "rotationPeriod", "rotPeriod", src=b)
            )
            age = try_float(get_first("age", "age_Myr", "age_million_years", src=b))
            mass = try_float(get_first("solarMasses", "mass", "stellarMass", src=b))
            temperature = try_float(
                get_first("surfaceTemperature", "temperature", "temp", src=b)
            )
            abs_magnitude = try_float(get_first("absoluteMagnitude", "absMag", src=b))

            if sr is None or rp is None or age is None:
                continue
            if isnan(sr) or isnan(rp) or isnan(age):
                continue
            if mass is not None and isnan(mass):
                mass = None
            if temperature is not None and isnan(temperature):
                temperature = None
            if abs_magnitude is not None and isnan(abs_magnitude):
                abs_magnitude = None

            # Determine spin (same logic as before)
            spin = "none"
            neutron_body_id = b.get("bodyId")

            if b.get("rings") or b.get("hasRings"):
                spin = "up"

            if spin == "none":
                for body in bodies:
                    body_parents = body.get("parents") or []
                    body_type = body.get("type") or ""
                    if body_type == "Planet" and body_parents:
                        first_parent = body_parents[0]
                        if (
                            "Star" in first_parent
                            and first_parent["Star"] == neutron_body_id
                        ):
                            spin = "up"
                            break

            if spin == "none":
                parents = b.get("parents") or []
                if parents and len(parents) > 0:
                    immediate_parent = parents[0]
                    if "Null" in immediate_parent:
                        barycenter_id = immediate_parent["Null"]
                        neutron_sma = try_float(
                            get_first("semiMajorAxis", "orbitalRadius", src=b)
                        )

                        for companion in bodies:
                            companion_parents = companion.get("parents") or []
                            if not companion_parents:
                                continue
                            comp_immediate_parent = companion_parents[0]

                            if "Null" in comp_immediate_parent:
                                comp_barycenter_id = comp_immediate_parent["Null"]
                                comp_body_id = companion.get("bodyId")

                                if (
                                    comp_barycenter_id == barycenter_id
                                    and comp_body_id != neutron_body_id
                                ):
                                    comp_type = companion.get("type") or ""

                                    if comp_type == "Star":
                                        comp_sma = try_float(
                                            get_first(
                                                "semiMajorAxis",
                                                "orbitalRadius",
                                                src=companion,
                                            )
                                        )
                                        comp_mass = try_float(
                                            get_first(
                                                "solarMasses",
                                                "mass",
                                                "stellarMass",
                                                src=companion,
                                            )
                                        )

                                        if (
                                            neutron_sma is not None
                                            and comp_sma is not None
                                        ):
                                            separation_ls = neutron_sma + comp_sma
                                            separation_sr = separation_ls * 215.03

                                            if separation_sr < 500:
                                                if (
                                                    mass is not None
                                                    and comp_mass is not None
                                                ):
                                                    if comp_mass > mass:
                                                        spin = "down"
                                                    else:
                                                        spin = "up"
                                                    break
                                                elif (
                                                    comp_mass is not None
                                                    and comp_mass > 0.5
                                                ):
                                                    if comp_mass > 1.4:
                                                        spin = "down"
                                                    else:
                                                        spin = "up"
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
                    "absoluteMagnitude": abs_magnitude,
                    "x": sys_x,
                    "y": sys_y,
                    "z": sys_z,
                    "spin": spin,
                }
            )
            stars_found += 1
            if max_points and len(points) >= max_points:
                return points

    print(f"Extracted {stars_found} neutron stars from systems")
    return points


def write_csv_by_region(points, output_dir="."):
    """Write neutron star data to region-specific CSV files"""
    start = time.time()
    print(f"Organizing {len(points)} points by region...")

    # Group points by region
    regions_dict = defaultdict(list)
    no_region_count = 0

    for point in points:
        x, y, z = point.get("x"), point.get("y"), point.get("z")
        if x is None or y is None or z is None:
            no_region_count += 1
            continue

        region_result = findRegion(x, y, z)
        if region_result is None:
            regions_dict["Unknown"].append(point)
        else:
            region_id, region_name = region_result
            regions_dict[f"Region_{region_id:02d}"].append(point)

    print(f"Found {len(regions_dict)} regions")
    if no_region_count > 0:
        print(f"Warning: {no_region_count} points had no coordinates")

    # Write each region to its own CSV
    fieldnames = [
        "system",
        "body",
        "solarRadius",
        "rotationalPeriod",
        "age",
        "mass",
        "temperature",
        "absoluteMagnitude",
        "x",
        "y",
        "z",
        "spin",
    ]

    files_written = []
    for region_name, region_points in sorted(regions_dict.items()):
        csv_path = os.path.join(output_dir, f"{region_name}.csv.gz")
        print(f"Writing {len(region_points)} points to {csv_path}")

        with gzip.open(csv_path, "wt", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for point in region_points:
                writer.writerow(
                    {
                        "system": point["system"],
                        "body": point["body"],
                        "solarRadius": point["solarRadius"],
                        "rotationalPeriod": point["rotationalPeriod"],
                        "age": point["age"],
                        "mass": point["mass"] if point["mass"] is not None else "",
                        "temperature": (
                            point["temperature"]
                            if point["temperature"] is not None
                            else ""
                        ),
                        "absoluteMagnitude": (
                            point["absoluteMagnitude"]
                            if point["absoluteMagnitude"] is not None
                            else ""
                        ),
                        "x": point["x"] if point["x"] is not None else "",
                        "y": point["y"] if point["y"] is not None else "",
                        "z": point["z"] if point["z"] is not None else "",
                        "spin": point["spin"],
                    }
                )
        files_written.append(csv_path)

    elapsed = time.time() - start
    print(f"Wrote {len(files_written)} region CSV files in {elapsed:.1f}s")
    return files_written


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
            "spin",
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
                    "spin": point["spin"],
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
                    "spin": row["spin"],
                }
            )

    elapsed = time.time() - start
    print(f"Read {len(points)} points from CSV in {elapsed:.1f}s")
    return points


def main():
    signal.signal(signal.SIGINT, signal_handler)

    p = argparse.ArgumentParser(
        description="Generate neutron star data CSV cache by region"
    )
    p.add_argument(
        "--input",
        "-i",
        default="/home/meddler/spansh/galaxy.json.gz",
        help="path to galaxy.json.gz",
    )
    p.add_argument(
        "--systems-jsonl",
        default="systems.jsonl",
        help="path to systems.jsonl cache file",
    )
    p.add_argument(
        "--output-dir",
        "-o",
        default=".",
        help="output directory for region CSV files",
    )
    p.add_argument("--max", type=int, default=0, help="max points to read (0 = all)")
    p.add_argument(
        "--progress-interval",
        type=int,
        default=200000,
        help="lines between progress messages (0 to disable)",
    )
    p.add_argument(
        "--force-reload",
        action="store_true",
        help="force reload from galaxy file even if systems.jsonl exists",
    )
    p.add_argument(
        "--track-region-ages",
        action="store_true",
        help="track average age of all stars by region and save to region_ages.json",
    )
    p.add_argument(
        "--region-ages-output",
        default="region_ages.json",
        help="output path for region ages JSON file",
    )

    args = p.parse_args()

    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    systems = None
    pts = None

    # Check if systems.jsonl exists
    if os.path.exists(args.systems_jsonl) and not args.force_reload:
        print(f"Found existing {args.systems_jsonl}, loading from cache...")
        systems = load_systems_jsonl(args.systems_jsonl)
        pts = extract_neutron_stars_from_systems(
            systems, max_points=(args.max if args.max > 0 else None)
        )
    else:
        print("Extracting neutron stars from galaxy file...")
        # extract_neutron_stars now writes systems.jsonl directly while streaming
        pts = extract_neutron_stars(
            args.input,
            max_points=(args.max if args.max > 0 else None),
            progress_interval=args.progress_interval,
            systems_jsonl_path=args.systems_jsonl,
            track_region_ages=args.track_region_ages,
            region_ages_path=args.region_ages_output,
        )

    if not pts:
        print("No neutron star data found.")
        return

    point_count = len(pts)
    print(f"Extracted {point_count} neutron stars")

    # Write region-specific CSV files
    print(f"Writing region CSV files to {args.output_dir}...")
    files = write_csv_by_region(pts, args.output_dir)

    print(
        f"Done: {point_count} neutron stars written to {len(files)} region files in {args.output_dir}"
    )
    print(f"Systems cache saved to {args.systems_jsonl} for faster future runs")


if __name__ == "__main__":
    main()
