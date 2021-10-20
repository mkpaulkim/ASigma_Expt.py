import numpy as np
import os
import hardware as hw
import pycubelib.general_functions as gf
import pycubelib.tkparts as tkp
import pycubelib.plotting_functions as pf
import pycubelib.files_functions as ff

pi = np.pi
pi2 = pi * 2
pi2limit = (-pi, pi)
cwd = os.path.dirname(__file__)
txtpath0 = 'E:/{{SeaGate}}/Dropbox/[[ PROJECTS.dbox ]]/project folders 2021/proj 2021-10 AlphaSigma/temp_data/aaa.txt'

# fp = tkp.tkwindow('AlphaSigma', (20, 20, 1100, 550))
fp = tkp.tkwindow('AlphaSigma', (20, 20, 1100, 350))
arduino = hw.Arduino()
qcam = hw.Quantalux()
pzt = hw.MDT694B()
ddr = hw.DDR25()
# print(f'> ASigma_expt: ready ...')

ent_tcam = tkp.ParamEntry(fp, (100, 30, 10), 300, 'T_cam, ms')
ent_texp = tkp.ParamEntry(fp, (100, 60, 10), 3.00, 't_exp, ms')
btn_cset = tkp.CmdButton(fp, (100, 90, 10), 'cam set', 'orange')
btn_cam = tkp.CmdButton(fp, (100, 120, 10), 'camera', 'dark orange')
ent_level = tkp.ParamEntry(fp, (100, 150, 10), 0, 'level', rw='r')
prog_level = tkp.ProgressBar(fp, (50, 180, 125), '')
ent_lam0 = tkp.ParamEntry(fp, (100, 240, 10), 0.6328, 'lam0, um')
ent_ax = tkp.ParamEntry(fp, (100, 270, 10), 10, 'ax, mm')
ent_roi = tkp.ParamEntry(fp, (50, 300, 20), f'{qcam.nx//2}, {qcam.ny//2}, 10, 10', 'roi')

ent_vpz = tkp.ParamEntry(fp, (300, 30, 10), 0, 'V pz')
btn_vset = tkp.CmdButton(fp, (300, 60, 10), 'set V pz', 'orange')
ent_v2pi = tkp.ParamEntry(fp, (300, 90, 10), 5.30, 'V 2pi')
btn_calib = tkp.CmdButton(fp, (300, 120, 10), 'calib V2pi', 'orange')
btn_getvv = tkp.CmdButton(fp, (300, 150, 10), 'get vv', 'orange')
ent_v1 = tkp.ParamEntry(fp, (300, 180, 10), 20.00, 'V1')
ent_nv = tkp.ParamEntry(fp, (300, 210, 10), 10, 'nV')
ent_n2pi = tkp.ParamEntry(fp, (300, 240, 10), 2, 'N2pi')
ent_v2 = tkp.ParamEntry(fp, (300, 270, 10), 0, 'V2', rw='r')
btn_hh = tkp.CmdButton(fp, (300, 300, 10), 'acq HH', 'dark orange')

ent_qddr = tkp.ParamEntry(fp, (500, 30, 10), 237.6000, 'Q_ddr')
btn_qset = tkp.CmdButton(fp, (500, 60, 10), 'set Q', 'orange')
ent_psi = tkp.ParamEntry(fp, (500, 90, 10), 18.0, 'psi')
ent_q0 = tkp.ParamEntry(fp, (500, 120, 10), 237.5959, 'Q0')
ent_q1 = tkp.ParamEntry(fp, (500, 150, 10), 237.6000, 'Q1')
btn_qns = tkp.CmdButton(fp, (500, 180, 10), 'get Qns', 'orange')
ent_lam12 = tkp.ParamEntry(fp, (500, 210, 10), 2000.0, 'lam_12, um')
ent_gamma = tkp.ParamEntry(fp, (500, 240, 10), 0.5, 'gamma')
ent_nh = tkp.ParamEntry(fp, (500, 270, 10), 5, 'Nh')
ent_lam1n_last = tkp.ParamEntry(fp, (500, 300, 10), 0, 'lam1n last', rw='r')

btn_madh = tkp.CmdButton(fp, (700, 30, 10), 'MADH', 'dark orange')
ent_qscan = tkp.ParamEntry(fp, (700, 90, 10), 0, 'Q scan', rw='r')
prog_qscan = tkp.ProgressBar(fp, (800, 90, 250), '')
ent_vscan = tkp.ParamEntry(fp, (700, 120, 10), 0, 'V scan', rw='r')
prog_vscan = tkp.ProgressBar(fp, (800, 120, 250), '')
btn_txt = tkp.CmdButton(fp, (700, 180, 10), 'txt path', 'orange')
ent_txtpath = tkp.ParamEntry(fp, (800, 180, 35), txtpath0, '')
ent_phspath = tkp.ParamEntry(fp, (800, 210, 35), '', 'phs_path', rw='r')
ent_txtpath.right()
ent_phspath.right()
# txtbox = tkp.TextBox(fp, (100, 350, 900, 150), '')
btn_adios = tkp.CmdButton(fp, (900, 30, 10), 'adios', 'indian red')


def main_loop():
    global roi

    if btn_cam.is_on():
        roi = tuple(np.float(a) for a in ent_roi.get_val(str).split(','))
        cc = show_cam()

    fp.after(tloop, main_loop)


def set_camera():
    global roi

    t_cam = ent_tcam.get_val()
    t_exp = ent_texp.get_val(typ=float)
    arduino.set_Tcam(t_cam)
    qcam.setup(t_exp=t_exp)
    roi = tuple(np.float(a) for a in ent_roi.get_val(str).split(','))

    # print(f'> set_cam: T_cam = {t_cam}, t_exp = {t_exp}')
    make_notes()


def show_cam():
    qcam.acquire_cc()
    capA = f'qcam.cc: iframe = {qcam.iframe}'
    cc = qcam.cc * 1.0
    pf.plotAAB(cc, figname='qcam.cc', capA=capA, ulimit=(0, 2**16-1), sxy=(.5, .5), aby=(3, 1), roi=roi, pause=.01)
    level = 100 * np.sum(cc) / (np.prod(np.shape(cc)) * (65536.0))
    ent_level.set_entry(f'{level:.1f}')
    prog_level.setval(level)


def set_vpz():
    vpz = ent_vpz.get_val(float)
    pzt.set_vpz(vpz)


def calib_v2pi():
    btn_calib.on()
    ss = []
    pzt.set_vpz(0)
    for iv, v in enumerate(vv):
        pzt.set_vpz(v)
        ent_vscan.set_entry(f'{v:.1f}')
        prog_vscan.setval(100 * (iv+1) / len(vv))
        qcam.acquire_cc()
        rx0, ry0, rx, ry = roi
        rr = qcam.cc[int(rx0-rx//2):int(rx0+rx//2), int(ry0-ry//2):int(ry0+ry//2)]
        ss += [np.sum(rr * 1.0) / np.prod(np.shape(rr))]
    v1, dv = (vv[0], vv[1] - vv[0])
    pf.graphB(ss, caption='ss vs vv', xpars=(v1, dv), sxy=(30, .3), line='-o')
    btn_calib.off()


def get_vv():
    global vv, phs

    v2pi = ent_v2pi.get_val(float)
    v1 = ent_v1.get_val(float)
    nv = ent_nv.get_val()
    n2pi = ent_n2pi.get_val()
    v2 = v1 + v2pi * n2pi
    ent_v2.set_entry(v2)
    dv = v2pi / nv
    vv = v1 + np.arange(nv * n2pi) * dv
    phs = pi2 * np.mod(np.arange(nv * n2pi), nv) / nv

    # print(f'> get_vv: V1 = {v1:.2f}, V_2pi = {v2pi}, NV = {nv}, N_2pi = {n2pi}, V2 = {v2:.2f}, dV = {dv:.3f}')
    print(gf.prn_list('vv', vv))
    print(gf.prn_list('phs', phs))
    make_notes()


def acq_hh(iq=1):
    btn_hh.on()
    qcam.acquire_cc()
    hh = np.array(qcam.cc) * (1j * 0.0)
    nv = len(vv)
    pzt.set_vpz(0)
    for iv, (v, ph) in enumerate(zip(vv, phs)):
        pzt.set_vpz(v, prn=False)
        ent_vscan.set_entry(f'{v:.1f}')
        prog_vscan.setval(100 * (iv+1) / len(vv))
        qcam.acquire_cc()
        hh += np.array(qcam.cc) * np.exp(1j * ph)

    # lam0 = ent_lam0.get_val(float)
    # k0 = pi2 / lam0
    # theta_n = 2 * (q_ns[iq] - q_ns[1]) * pi/180
    # gg = np.exp(1j * k0 * np.sin(theta_n) * xx)
    # hh = hh * gg

    hha = np.abs(hh) / nv / (2 ** 14)
    hhp = np.angle(hh)
    pf.plotAAB(hha, 'HH_abs', 'HH_abs', roi=roi, sxy=(.3, .3), pause=0.1)
    pf.plotAAB(hhp, 'HH_phs', 'HH_phs', roi=roi, sxy=(.3, .3), pause=0.1)
    btn_hh.off()
    return hha, hhp


def set_ddr_q():
    q_ddr = ent_qddr.get_val(float)
    ddr.goto(q_ddr)


def get_qns():
    global q_ns, lam_ns, lam_1ns
    d2r = pi/180
    r2d = 180/pi

    psi = ent_psi.get_val(float)
    q0 = ent_q0.get_val(float)
    q1 = ent_q1.get_val(float)
    lam0 = ent_lam0.get_val(float)
    lam12 = ent_lam12.get_val(float)
    gamma = ent_gamma.get_val(float)
    nh = ent_nh.get_val()
    theta1 = psi + 2. * (q1 - q0)
    lam1 = lam0 * np.cos(theta1 * d2r)
    lam_ns = [0, lam1]
    lam_1ns = [0, 0]
    q_ns = [0, q1]
    for n in range(2, nh+1):
        lam1n = gamma**(n-2) * lam12
        lamn = 1./(1./lam1n + 1./lam1)
        thetan = np.arccos(lamn / lam0) * r2d
        qn = q0 + (thetan - psi) / 2.
        lam_ns += [lamn]
        lam_1ns += [lam1n]
        q_ns += [qn]
    ent_lam1n_last.set_entry(f'{lam_1ns[-1]:.1f}')

    # gf.prn_list('lam_1ns', lam_1ns, 1)
    # gf.prn_list('q_ns', q_ns, 4)
    # gf.prn_list('lam_ns', lam_ns, 8)
    make_notes()


def madh():
    btn_madh.on()
    # abslimit = (0, 2**14-1)
    lam0 = ent_lam0.get_val(float)
    k0 = pi2/lam0

    ff.write_txt(make_notes(), ent_txtpath.get_val(str))

    shha = np.copy(blank)
    for iq, q in enumerate(q_ns[1:]):
        ddr.goto(q)
        ent_qscan.set_entry(f'{q:.4f}')
        prog_qscan.setval(100 * (iq+1) / len(q_ns[1:]))
        hha, hhp = acq_hh(iq+1)

        theta_n = 2 * (q_ns[iq+1] - q_ns[1]) * pi/180
        ggp = k0 * np.sin(theta_n) * xx
        hhp = np.mod(hhp + ggp + pi, pi2) - pi

        phs_path = ent_txtpath.get_val(str).replace('.txt', f'_{iq+1}p.png')
        ent_phspath.set_entry(phs_path)
        ff.write_png(hhp, phs_path, pi2limit)
        shha += hha

    abs_path = ent_txtpath.get_val(str).replace('.txt', f'_aa.png')
    hha = shha / len(q_ns[1:])
    ent_phspath.set_entry(abs_path)
    ff.write_png(hha, abs_path, (0., 1.))
    btn_madh.off()


def set_txtpath():
    from tkinter import filedialog

    txt_path = filedialog.asksaveasfilename(title='HH.txt file path', filetypes=[('txt files', '*.txt')])
    if txt_path.find('.txt') < 0:
        txt_path += '.txt'
    ent_txtpath.set_entry(txt_path)

    make_notes()


def make_notes():
    cam_notes = f'> Tcam_ms = {ent_tcam.get_val()}; ' \
                f't_exp_ms = {ent_texp.get_val(float)}; ' \
                f'lam0_um = {ent_lam0.get_val(float)}; ' \
                f'ax_mm = {ent_ax.get_val(float)}; \n'
    vpz_notes = f'> V_2pi = {ent_v2pi.get_val(float)}; ' \
                f'V1 = {ent_v1.get_val(float)}; ' \
                f'nV = {ent_nv.get_val()}; ' \
                f'N2pi = {ent_n2pi.get_val()} ' \
                f'V2 = {ent_v2.get_val(float)}; \n'
    qq_notes = f'> psi = {ent_psi.get_val(float)}; ' \
               f'Q0 = {ent_q0.get_val(float)}; ' \
               f'Q1 = {ent_q1.get_val(float)}; ' \
               f'lam_12_um = {ent_lam12.get_val(float)}; ' \
               f'gamma = {ent_gamma.get_val(float)}; ' \
               f'Nh = {ent_nh.get_val()}; ' \
               f'lam_1N = {ent_lam1n_last.get_val(float)}; \n'
    list_notes = gf.prn_list('q_ns', q_ns, 4) + '\n' + \
                 gf.prn_list('lam_ns', lam_ns, 8) + '\n' + \
                 gf.prn_list('lam_1ns', lam_1ns, 1) + '\n'
    header0 = f'> [{gf.path_parts(ent_txtpath.get_val(str), 3)[0]}] \n'
    header1 = f'> >> nx = {qcam.nx}; ny = {qcam.ny}; dx_um = {ent_ax.get_val(float) / qcam.nx:.3f}; ' \
             f'nh = {ent_nh.get_val()}; '+gf.prn_list('lam_ns', lam_ns)[2:] + '; \n'
    notes = header0 + cam_notes + vpz_notes + qq_notes + list_notes + header1

    print(f'\n{notes}')
    # txtbox.update(notes)
    return notes


def adios():
    qcam.dispose()
    fp.destroy()
    print(f'> adios amigos ...')
    quit()


btn_cset.command(set_camera)
btn_cam.command(btn_cam.switch)
btn_vset.command(set_vpz)
btn_calib.command(calib_v2pi)
btn_getvv.command(get_vv)
btn_hh.command(acq_hh)
btn_qset.command(set_ddr_q)
btn_qns.command(get_qns)
btn_madh.command(madh)
btn_txt.command(set_txtpath)
btn_adios.command(adios)

vv = []
q_ns = []
lam_ns = []
lam_1ns = []

set_camera()
set_vpz()
get_vv()
set_ddr_q()
get_qns()
btn_cam.on()

blank = qcam.cc * 0.0
nx, ny = qcam.nxy
ax = ent_ax.get_val(float) * 1e3
dx = ax / nx
ay = dx * ny
x = np.arange(nx) * dx - ax/2
y = np.arange(ny) * dx - ay/2
xx, yy = np.meshgrid(x, y)

tloop = 10
fp.after(tloop, main_loop)
fp.mainloop()



