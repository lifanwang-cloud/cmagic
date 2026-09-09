"""The standardized diagnostic panel: color-magnitude locus, window boundaries,
fit and mode, color-curve inset, annotation block."""
import numpy as np


def make_panel(P, path, title=''):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 6))
    ok = P.get('ok', False) and P.get('failed') is None
    loc = P.get('locus')
    if loc and loc['phases']:
        ph = np.array(loc['phases']); col = np.array(loc['colors'])
        B = np.array(loc['B'])
        ax.plot(col, B, 'o', color='grey', ms=5, alpha=0.55, label='paired nights')
        s = P.get('sel')
        if s and s['phases']:
            sc = np.array(s['colors']); sB = np.array(s['B'])
            ax.plot(sc, sB, 'o', color='tab:blue', ms=7, mec='k',
                    label=f"selected (mode {P.get('mode', '?')})")
            for p_, c_, b_ in zip(s['phases'], sc, sB):
                ax.annotate(f'+{p_:.0f}', (c_, b_), textcoords='offset points',
                            xytext=(5, 5), fontsize=7, color='tab:blue')
            xx = np.linspace(min(col.min(), 0.1), col.max() + 0.1, 40)
            if P.get('mode') == 'Q' and 'quad_coefs' in P:
                a0, a1 = P['quad_coefs'][0], P['quad_coefs'][1]
                a2 = P['quad_coefs'][2] if len(P['quad_coefs']) > 2 else 0.
                ax.plot(xx, a0 + a1 * (xx - 0.6) + a2 * (xx - 0.6) ** 2, 'r-',
                        lw=1.5, label=f"quadratic, rms {P.get('rms')}")
                ax.plot(0.6, a0, '*', color='gold', ms=17, mec='k')
            elif P.get('mode') == 'R' and P.get('m_star') is not None:
                cs = P['c_star']; bu = P.get('beta_used', P.get('beta_fixed', 1.94))
                ax.plot(xx, P['m_star'] + bu * (xx - cs), 'r-', lw=1.5,
                        label=f'beta={bu} at c*={cs}')
                ax.plot(cs, P['m_star'], '*', color='gold', ms=17, mec='k',
                        label=f"m*={P['m_star']}")
            elif P.get('B_BV06_raw') is not None:
                bf = P.get('beta_fixed', 1.94)
                ax.plot(xx, P['B_BV06_raw'] + bf * (xx - 0.6), 'r-', lw=1.5,
                        label=f"beta={bf} forced, rms {P.get('rms')}")
                ax.plot(0.6, P['B_BV06_raw'], '*', color='gold', ms=17, mec='k')
            if P.get('beta_free_valid'):
                b0 = P.get('B_BV06_raw', P.get('m_star'))
                x0 = 0.6 if P.get('B_BV06_raw') is not None else P.get('c_star', 0.6)
                ax.plot(xx, b0 + P['beta_free'] * (xx - x0), ':', color='purple',
                        lw=1.4,
                        label=f"free beta={P['beta_free']}+-{P['ebeta_free']}")
            ax.axvline(sc.min(), ls='--', lw=1.1, color='green',
                       label=f"entry (window +{P.get('window', ['?'])[0]} d)")
            ax.axvline(sc.max(), ls='--', lw=1.1, color='firebrick',
                       label=f"exit (+{P.get('window', ['?', '?'])[1]} d)")
        ax.invert_yaxis()
        ins = ax.inset_axes([0.63, 0.08, 0.34, 0.28])
        o_ = np.argsort(ph)
        ins.plot(ph[o_], col[o_], '.-', color='grey', ms=3, lw=0.7)
        if s and s['phases']:
            ins.plot(s['phases'], s['colors'], 'o', color='tab:blue', ms=3)
        if P.get('t_BVmax') is not None:
            ins.axvline(P['t_BVmax'], color='firebrick', ls='--', lw=0.9)
        if P.get('window'):
            ins.axvspan(P['window'][0], P['window'][1], color='gold', alpha=0.2)
        ins.set_xlabel('phase [d]', fontsize=6)
        ins.set_ylabel('B-V', fontsize=6)
        ins.tick_params(labelsize=6)
        ins.set_title('color curve + window', fontsize=6.5)
    lines = [f"source: {P.get('source')}",
             f"mode {P.get('mode', '?')}; window {P.get('window')} d"]
    if P.get('n_sel'):
        lines.append(f"n={P['n_sel']} ph {P.get('ph_span')} col {P.get('col_span')}")
    for k in ('B_BV06_raw', 'm_star', 'c_star', 'E_guess', 'E_host', 'B_BV_final',
              'M_BV', 'mu', 'D_mpc'):
        if P.get(k) is not None:
            lines.append(f'{k}={P[k]}')
    if P.get('beta_free') is not None:
        lines.append(f"beta_free={P['beta_free']}+-{P.get('ebeta_free')} "
                     f"[{'valid' if P.get('beta_free_valid') else 'diagnostic only'}]")
    for f in P.get('flags', []):
        lines.append(str(f)[:60])
    if not ok:
        lines.append('FAILED: ' + str(P.get('failed')))
    ax.text(0.02, 0.02, '\n'.join(lines), transform=ax.transAxes, fontsize=7.2,
            va='bottom', bbox=dict(facecolor='ivory', edgecolor='grey', alpha=0.9))
    ax.set_xlabel('B - V [mag]'); ax.set_ylabel('B [mag]')
    ax.set_title(title or 'CMAGIC diagnostic', fontsize=11,
                 color='k' if ok else 'tab:red')
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc='upper right')
    fig.tight_layout(); fig.savefig(path, dpi=170); plt.close(fig)
    return path
