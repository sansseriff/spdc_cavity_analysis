import * as math from "mathjs";


export function linspace(start, end, num) {
    return math
        .range(start, end, (end - start) / (num - 1), true)
        .toArray();
}

export function Sellimeier_PPLN(WL, pol, T) {
    WL = WL.map((wl) => wl / 1e-6); // in the Paper edwards 1984 Wl is menitoned in um

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