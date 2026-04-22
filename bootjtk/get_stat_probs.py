from scipy.stats import kendalltau as kt
from scipy.stats import circmean as sscircmean
from scipy.stats import circstd as sscircstd
import numpy as np


def get_stat_probs(dorder, new_header, triples, dref, size):
    d_taugene, d_pergene, d_phgene, d_nagene = {}, {}, {}, {}
    rs = []
    for kkey in dorder:
        res = np.array([get_matches(kkey, triple, dref, new_header) for triple in triples])
        r = pick_best_match(res)
        d_taugene[r[0]] = d_taugene.get(r[0], 0) + dorder[kkey]
        d_pergene[r[2]] = d_pergene.get(r[2], 0) + dorder[kkey]
        d_phgene[r[3]] = d_phgene.get(r[3], 0) + dorder[kkey]
        d_nagene[r[4]] = d_nagene.get(r[4], 0) + dorder[kkey]
        count = int(np.round(size * dorder[kkey]))
        rs.extend([r] * count)
    rs = np.array(rs)
    m_tau = np.mean(rs[:, 0])
    s_tau = np.std(rs[:, 0])
    m_per = np.mean(rs[:, 2])
    s_per = np.std(rs[:, 2])
    m_ph = sscircmean(rs[:, 3], high=24, low=0)
    s_ph = sscircstd(rs[:, 3], high=24, low=0)
    m_na = sscircmean(rs[:, 4], high=24, low=0)
    s_na = sscircstd(rs[:, 4], high=24, low=0)
    return [m_per, s_per, m_ph, s_ph, m_na, s_na], [m_tau, s_tau], d_taugene, d_pergene, d_phgene, d_nagene


def generate_base_reference(header, waveform="cosine", period=24., phase=0., width=12.):
    ZTs = np.array(header, dtype=float)
    coef = 2.0 * np.pi / period
    w = (width * coef) % (2.0 * np.pi)
    tpoints = ((ZTs - phase) * coef) % (2.0 * np.pi)
    if waveform == 'cosine':
        return np.where(
            tpoints <= w,
            np.cos(tpoints / (w / np.pi)),
            np.cos((tpoints + 2.0 * (np.pi - w)) * np.pi / (2.0 * np.pi - w))
        )
    elif waveform == 'trough':
        return np.where(
            tpoints <= w,
            1.0 - tpoints / w,
            (tpoints - w) / (2.0 * np.pi - w)
        )
    elif waveform == 'impulse':
        d = np.minimum(tpoints, np.abs(2.0 * np.pi - tpoints))
        return np.maximum(-2.0 * d / (3.0 * np.pi / 4.0) + 1.0, 0.0)
    elif waveform == 'step':
        return np.where(tpoints < np.pi, 1.0, 0.0)


def farctanh(x):
    return float(np.arctanh(np.clip(x, -0.99, 0.99)))


def periodic(x):
    x = float(x)
    while x > 12:
        x -= 24.
    while x <= -12:
        x += 24.
    return x


def pick_best_match(res):
    res = np.array(res)
    taus = res[:, 0]
    tau_mask = (max(taus) == taus)
    if np.sum(tau_mask) == 1:
        return res[int(np.argmax(tau_mask))]

    res = res[tau_mask]
    phases = np.abs(res[:, 3] - res[:, 5])
    phasemask = (min(phases) == phases)
    if np.sum(phasemask) == 1:
        return res[int(np.argmax(phasemask))]

    res = res[phasemask]
    diffs = np.abs(res[:, 4] - res[:, 6])
    diffmask = (min(diffs) == diffs)
    if np.sum(diffmask) == 1:
        return res[int(np.argmax(diffmask))]

    return res[np.random.randint(len(res))]


def get_waveform_list(periods, phases, widths):
    triples = []
    for period in periods:
        seen = set()
        for phase in phases:
            for width in widths:
                nadir = (phase + width) % period
                # deduplicate waveforms where peak and nadir positions are swapped
                key = (float(nadir), float(phase))
                if key not in seen:
                    seen.add(key)
                    triples.append([float(period), float(phase), float(width)])
    return np.array(triples, dtype=float)


def make_references(new_header, triples, waveform='cosine'):
    dref = {}
    for triple in triples:
        period, phase, width = triple
        reference = generate_base_reference(new_header, waveform, period, phase, width)
        dref[(period, phase, width)] = reference
    return dref


def get_matches(kkey, triple, d_ref, new_header):
    reference = d_ref[tuple(triple)]
    period, phase, width = triple
    nadir = (phase + width) % period
    tau, p = kt(reference, kkey)
    p = p / 2.0
    tau = farctanh(tau)
    maxloc = new_header[kkey.index(max(kkey))]
    minloc = new_header[kkey.index(min(kkey))]
    r = [tau, p, period, phase, nadir, maxloc, minloc]
    if tau < 0:
        r = [abs(tau), p, period, nadir, phase, maxloc, minloc]
    return r
