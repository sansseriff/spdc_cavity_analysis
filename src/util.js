import { create } from "d3";
import * as math from "mathjs";


export function linspace(start, end, num) {
    return math
        .range(start, end, (end - start) / (num - 1), true)
        .toArray();
}

export function arange(start, stop, step) {
    if (typeof stop === 'undefined') {
        // one param defined
        stop = start;
        start = 0;
    }

    if (typeof step === 'undefined') {
        step = 1;
    }

    if ((step > 0 && start >= stop) || (step < 0 && start <= stop)) {
        return [];
    }

    const result = [];
    for (let i = start; step > 0 ? i < stop : i > stop; i += step) {
        result.push(i);
    }

    return result;
}

export function Sellimeier_PPLN(WL, pol, T) {
    WL = WL.map((wl) => wl * 1e6); // in the Paper edwards 1984 Wl is menitoned in um

    const T0 = 24.5;
    const F = (T - T0) * (T + T0 + 546);

    let A1, A2, A3, A4, B1, B2, B3;

    if (pol === "eray") {
        A1 = 4.9048;
        A2 = 0.11775;
        A3 = 0.21802;
        A4 = 0.027153;
        B1 = 2.2314e-8;
        B2 = -2.9671e-8;
        B3 = 2.1429e-8;
    } else if (pol === "oray") {
        A1 = 4.582;
        A2 = 0.09921;
        A3 = 0.2109;
        A4 = 0.02194;
        B1 = 5.2716e-8;
        B2 = -4.91431e-8;
        B3 = 2.2971e-7;
    }

    const n = WL.map((wl) =>
        Math.sqrt(
            A1 +
                (A2 + B1 * F) /
                    (Math.pow(wl, 2) - Math.pow(A3 + B2 * F, 2)) +
                B3 * F -
                A4 * Math.pow(wl, 2),
        ),
    );
    return n;
}

export function sellimeier_ppln_single(WL, pol, T) {
    WL = WL * 1e6; // in the Paper edwards 1984 Wl is menitoned in um

    const T0 = 24.5;
    const F = (T - T0) * (T + T0 + 546);

    let A1, A2, A3, A4, B1, B2, B3;

    // console.log("F: ", F)


    if (pol === "eray") {
        A1 = 4.9048;
        A2 = 0.11775;
        A3 = 0.21802;
        A4 = 0.027153;
        B1 = 2.2314e-8;
        B2 = -2.9671e-8;
        B3 = 2.1429e-8;
    } else if (pol === "oray") {
        A1 = 4.582;
        A2 = 0.09921;
        A3 = 0.2109;
        A4 = 0.02194;
        B1 = 5.2716e-8;
        B2 = -4.91431e-8;
        B3 = 2.2971e-7;
    }

    const n = Math.sqrt(
        A1 +
            (A2 + B1 * F) /
                (Math.pow(WL, 2) - Math.pow(A3 + B2 * F, 2)) +
            B3 * F -
            A4 * Math.pow(WL, 2),
    );

    return n;
}

// export function Sellimeier_PPKTP(n, WL) {
//     WL = WL * 1e6;

//     let A, B, C, D;
//     if (n === "nx") {
//         A = 2.1146;
//         B = 0.89188;
//         C = 0.20861;
//         D = 0.01320;
//     } else if (n === "ny") {
//         A = 2.1518;
//         B = 0.87862;
//         C = 0.21801;
//         D = 0.01327;
//     } else if (n === "nz") {
//         A = 2.3136;
//         B = 1.00012;
//         C = 0.23831;
//         D = 0.01679;
//     }

//     n = Math.sqrt(A + B / (1 - Math.pow(C / WL, 2)) - D * Math.pow(WL, 2));

//     return n;
// }

export function Sellimeier_PPKTP(direction, WL) {
    // this is acctually the nx value, the curve matches the data sent by raicol

    // nx
    // return math.sqrt(3.29100 + (0.04140/(WL**2 - 0.03978)) + 9.35522/(WL**2 - 31.45571))


    // ny 
    // return math.sqrt(3.45018 + 0.04341/(WL**2 - 0.04597) + 16.98825/(WL**2 - 39.43799))


    // nz
    return math.sqrt(4.59423 + 0.06206/(WL**2 - 0.04763) + 110.80672/(WL**2 - 86.12171))
}


export function n_raicol_ppktp(WL, T, debug=false) {
    if (debug) {
        console.log("n_raicol_ppktp WL: ", WL, "n: ", Sellimeier_PPKTP("nz", WL));
    }

    return Sellimeier_PPKTP("nz", WL);
    // return 1.9
}

export function n_ppln(WL, T, debug=false) {
    if (debug) {
        console.log("n_ppln WL: ", WL, "n: ", sellimeier_ppln_single(WL, "eray", T));
    }
    return sellimeier_ppln_single(WL, "eray", T)
}

export function AiryFunction(R1, R2, L, neff, WL) {
    if (neff.length !== WL.length) {
        throw new Error('neff and WL must be the same length');
    }

    const A = neff.map((n, i) => {
        const wl = WL[i];
        const phi = (2 * Math.PI * L * n) / wl;
        return Math.pow(
            1 +
                (4 * Math.sqrt(R1 * R2) * Math.pow(Math.sin(phi), 2)) /
                    Math.pow(1 - Math.sqrt(R1 * R2), 2),
            -1,
        );
    });

    return A;
}

export function finesse(R1, R2, alpha, L) {
    const rho = R1 * R2 * Math.pow(10, (-2 * alpha * L) / 10);
    const Finesse_val =
        Math.PI /
        (2 * Math.asin((1 - Math.sqrt(rho)) / (2 * Math.pow(rho, 0.25))));
    return Finesse_val;
}

export function phase_matching(WL1, WL2, index_function, T, polling_period, length, debug_show = false) {
    const delK = delta_k(WL1, WL2, index_function, T, polling_period, debug_show);
    let phase = Math.pow((Math.sin(delK*length/2)/(delK*length/2)),2);

    // if (debug_show) {
    //     phase = 20;
    // }
    return phase;
}



export function wl_pump_by_energy_conservation(WL1, WL2, index_function, T) {
    const nWL1 = index_function(WL1, T);
    const nWL2 = index_function(WL2, T);
    // const pumpWL = 1/((nWL1/WL1) + (nWL2/WL2));
    const pumpWL = 1/((1/WL1) + (1/WL2));
    return pumpWL;
}

export function delta_k(WL1, WL2, index_function, T, polling_period, debug_show = false) {
    // if (debug_show) {
    //     console.log("WL1: ", WL1);
    //     console.log("WL2: ", WL2);
    // }

    const nWL1 = index_function(WL1, T, debug_show);
    const nWL2 = index_function(WL2, T, debug_show);

    const WL_pump = wl_pump_by_energy_conservation(WL1, WL2, index_function, T);
    const nPump = index_function(WL_pump, T, debug_show);

    

    const delta_k = 2 * Math.PI * (nPump/WL_pump - nWL1/WL1 - nWL2/WL2 - 1/polling_period);

    if (debug_show) {
        // console.log(nPump)
        // console.log("nWL1: ", nWL1);
        // console.log("nWL2: ", nWL2);
        // console.log("nPump: ", nPump);
        // console.log("pump wavelength: ", WL_pump);
        // console.log(WL1, WL2, WL_pump);
        // // console.log("WL1", WL1);
        // console.log("nPump/WL_pump", nPump/WL_pump);
        // console.log("nWL1/WL1", nWL1/WL1);
        // console.log("nWL2/WL2", nWL2/WL2);
        // console.log("polling_period", polling_period);
        // console.log("1/polling_period", 1/polling_period);
        // console.log("nPump/WL_pump - nWL1/WL1 - nWL2/WL2", nPump/WL_pump - nWL1/WL1 - nWL2/WL2);
        // console.log(delta_k)
    }
    return delta_k
}


export function WL_by_energy_conservation(WL1, WL2) {
    const nWL1med = Sellimeier_PPLN(WL1, "eray", 24.5);
    const nWL2med = Sellimeier_PPLN(WL2, "eray", 24.5);

    const WL3med = WL1.map((wl1, i) =>
        Math.abs(1 / (nWL1med[i] / wl1 - nWL2med[i] / WL2[i])),
    );
    const WL3vacuum = WL1.map((wl1, i) =>
        Math.abs(1 / (1 / wl1 - 1 / WL2[i])),
    );
    const nWL3med = Sellimeier_PPLN(WL3vacuum, "eray", 24.5);

    return [WL3vacuum, WL3med];
}

export function create_50ghz_wl() {
    const freqs = arange(191.15e12, 196.15e12, 50e9);
    // console.log("frequencies: ", freqs);
    const c = 299792458;
    const wl = freqs.map((f) => c*1e9 / f);
    // console.log("wavelengths: ", wl);
    const midPoint = Math.ceil(wl.length / 2);
    // const wl1 = wl.slice(0, midPoint);
    // const wl2 = wl.slice(midPoint);
    
    // console.log("First half of wavelengths: ", wl1);
    // console.log("Second half of wavelengths: ", wl2);
    return wl.reverse();
}

export function interp(x, xp, yp) {
    if (xp.length !== yp.length) {
        throw new Error('xp and yp must have the same length');
    }

    if (xp.length === 1) {
        return Array.isArray(x) ? Array(x.length).fill(yp[0]) : yp[0];
    }

    const interpolate = (xi) => {
        for (let i = 0; i < xp.length - 1; i++) {
            if (xi >= xp[i] && xi <= xp[i + 1]) {
                const slope = (yp[i + 1] - yp[i]) / (xp[i + 1] - xp[i]);
                return yp[i] + slope * (xi - xp[i]);
            }
        }

        // // Slope Extrapolation
        // if (xi < xp[0]) {
        //     return yp[0] + (xi - xp[0]) * (yp[1] - yp[0]) / (xp[1] - xp[0]);
        // } else {
        //     const n = xp.length;
        //     return yp[n - 1] + (xi - xp[n - 1]) * (yp[n - 1] - yp[n - 2]) / (xp[n - 1] - xp[n - 2]);
        // }

        // Constant value extrapolation
        if (xi < xp[0]) {
            return yp[0];
        } else {
            return yp[yp.length - 1];
        }
    };

    return Array.isArray(x) ? x.map(interpolate) : interpolate(x);
}


export function RaicolNDataPPKTP(data) {
    this.wl = data.x.map((wl) => wl * 1e-6);
    this.n = data.y


    // I need to use an arrow function here because I end up not 
    // calling the function directly 'on' the object
    this.index_function = (WL, T, debug=false) => {
        // console.log("this.wl: ", this.wl)

        const delta_T = T - 24;

        // if (debug) {
        //     console.log("WL: ", WL)
        // }
        const wl = WL * 1e6;

        const del_n_del_t = ((0.1717/(wl**3)) - (0.5353/(wl**2)) + (0.8416/wl) + 0.1627) * 1e-5;
        // console.log(del_n_del_t)
        // if (debug) {
        //     console.log("del_n_del_t: ", del_n_del_t)
        // }

        return Number(Sellimeier_PPKTP("nx", WL*1e6)) + del_n_del_t * delta_T;
        // return interp(WL, this.wl, this.n);
    }
}

function closestIndex(sortedArray, target) {
    let start = 0;
    let end = sortedArray.length - 1;

    while (end - start > 1) {
        let mid = Math.floor((start + end) / 2);
        if (sortedArray[mid] < target) {
            start = mid;
        } else {
            end = mid;
        }
    }

    if (target - sortedArray[start] <= sortedArray[end] - target) {
        return start;
    }
    return end;
}

export function FilterComb(filter_spectrum, wl_array) {
    this.x = filter_spectrum.x;
    this.y = filter_spectrum.y;
    this.wl_array = wl_array;
    this.out_array = Array(wl_array.length).fill(0);
    console.log("length of out_array: ", this.out_array.length)

    this.spacing_wl = create_50ghz_wl();
    const center_of_filter_idx = 496;

    const range = 0.2

    for (let i = 0; i < this.spacing_wl.length; i++) {
        
        const spacing = this.spacing_wl[i];

        const current_filter_xp = this.x.map((wl) => wl + spacing); // this makes x[496] equal to spacing
        const current_filter_yp = this.y;

        // index of wl in 'wl_array' that is closest to 'spacing'
        const closest_wl_idx = closestIndex(this.wl_array, spacing);

        this.out_array[closest_wl_idx] = interp(wl_array[closest_wl_idx], current_filter_xp, current_filter_yp);

        for (let j = 1; j < this.wl_array.length; j++) {
            const raw_out_idx_left = closest_wl_idx - j;
            const raw_out_idx_right = closest_wl_idx + j;
            const max_idx = this.wl_array.length - 1;
            const min_idx = 0;



            if ((raw_out_idx_right > max_idx) || (raw_out_idx_left > max_idx)|| (raw_out_idx_left < 0 || (raw_out_idx_right < 0)) ) {
                break
            }
            this.out_array[raw_out_idx_left] = interp(wl_array[raw_out_idx_right], current_filter_xp, current_filter_yp);
            this.out_array[raw_out_idx_right] = interp(wl_array[raw_out_idx_right], current_filter_xp, current_filter_yp);

            

            if (((wl_array[closest_wl_idx] - wl_array[raw_out_idx_left]) > range) || ((wl_array[raw_out_idx_right] - wl_array[closest_wl_idx]) > range)) {
                break;
            }
        }
    }
    console.log("length of out_array: ", this.out_array.length)
}