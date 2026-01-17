// This file provides the jet prediction function for use in the HTML plot
// Copied from user-provided TypeScript, converted to JavaScript for browser use

/**
 * Predicts jet height (Ls) and half-angle (deg) for a body
 * @param {object|number} body - Body object with properties or massSolar if using individual params
 * @param {number} [radiusSolar] - Solar radii (if using individual params)
 * @param {number} [rotationalPeriodDays] - Days (if using individual params)
 * @returns {{ jetHeightLs: number, halfAngleDeg: number }}
 */
function predictJet(body, radiusSolar, rotationalPeriodDays) {
    // Handle both object and individual parameter formats
    let massSolar, radiusSolarVal, rotationalPeriodDaysVal;
    
    if (typeof body === 'object' && body !== null) {
        // Object format: predictJet({massSolar: 1.5, radiusSolar: 0.012, rotationalPeriodDays: 0.002})
        massSolar = body.massSolar;
        radiusSolarVal = body.radiusSolar;
        rotationalPeriodDaysVal = body.rotationalPeriodDays;
        console.log('[predictJet] Object format - Mass:', massSolar, 'Radius:', radiusSolarVal, 'Period:', rotationalPeriodDaysVal);
    } else {
        // Individual parameters format: predictJet(1.5, 0.012, 0.002)
        massSolar = body;
        radiusSolarVal = radiusSolar;
        rotationalPeriodDaysVal = rotationalPeriodDays;
        console.log('[predictJet] Individual params - Mass:', massSolar, 'Radius:', radiusSolarVal, 'Period:', rotationalPeriodDaysVal);
    }
    // Constants
    const G = 6.6743e-11; // m^3/kg/s^2
    const c = 3e8;        // m/s
    const M_sun = 1.98847e30;  // kg
    const R_sun = 6.957e8;     // m

    // Convert input to SI
    const M = massSolar * M_sun;
    const R = radiusSolarVal * R_sun;
    const T = rotationalPeriodDaysVal * 24 * 3600; // seconds
    
    console.log('[predictJet] SI units - M:', M, 'kg, R:', R, 'm, T:', T, 's');

    // Compactness
    const C = (G * M) / (R * c ** 2);

    // Surface gravity
    const g = (G * M) / (R ** 2);

    // Centrifugal acceleration at equator
    const aC = (4 * Math.PI ** 2 * R) / (T ** 2);

    // Effective gravity
    const gEff = g - aC;

    // Jet height [Ls] scaling constants
    const k1 = 20; // scales inverse compactness
    const k2 = 2;  // rotation boost factor

    let jetHeightLs = k1 * (1 / C) * (1 + k2 * (aC / g));
    jetHeightLs = Math.max(0.5, Math.min(jetHeightLs, 25));

    // Half-angle [deg] - Physics-based model using compactness and rotation
    const thetaMin = 3;
    const thetaMax = 45;
    
    // Compactness-based model: more compact objects have narrower jets
    // Typical neutron star compactness: 0.1-0.3
    // Model: theta = thetaMax * exp(-alpha * C) + thetaMin
    const alpha = 8.0; // tuning parameter for compactness sensitivity
    let compactnessAngle = thetaMax * Math.exp(-alpha * C) + thetaMin;
    
    // Rotation effect: faster rotation can focus jets (up to a limit)
    // Breakup frequency for neutron star: f_K = sqrt(GM/R^3) / (2*pi)
    const f_K = Math.sqrt(G * M / (R ** 3)) / (2 * Math.PI); // Hz
    const f_rot = 1 / T; // rotation frequency Hz
    const f_ratio = Math.min(f_rot / f_K, 0.8); // limit to 80% of breakup
    
    // Rotation focusing factor: 0.7-1.0 (faster rotation = more focused)
    const rotationFactor = 1.0 - 0.3 * f_ratio;
    
    let halfAngleDeg = compactnessAngle * rotationFactor;
    halfAngleDeg = Math.max(thetaMin, Math.min(halfAngleDeg, thetaMax));

    const result = { jetHeightLs, halfAngleDeg };
    console.log('[predictJet] Result - Jet Height:', jetHeightLs.toFixed(2), 'Ls, Half Angle:', halfAngleDeg.toFixed(1), '°');
    return result;
}

// Test function to verify the jet prediction works
function testJetPrediction() {
    console.log('[testJetPrediction] Testing jet prediction with example data...');
    
    const body = {
        massSolar: 1.5,          // 1.5 solar masses
        radiusSolar: 0.012,      // 0.012 solar radii (~neutron star)
        rotationalPeriodDays: 0.002 // 0.002 days (~3 minutes)
    };
    
    const result = predictJet(body);
    console.log(`[testJetPrediction] Predicted jet height: ${result.jetHeightLs.toFixed(2)} Ls`);
    console.log(`[testJetPrediction] Predicted half-angle: ${result.halfAngleDeg.toFixed(1)}°`);
    
    return result;
}
