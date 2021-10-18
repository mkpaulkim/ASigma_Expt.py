import numpy as np
import matplotlib.pyplot as plt

dpi = 100
font0 = 'Consolas 10'


def plotAA(AA, figname='plotAA', caption='', xypars=(0, 0, 1, 1), ulimit=(), sxy=(1, 1), cmap='gray', pause=0.):
    """
    xypars = (x0, y0, dx, dy); () no axes;
    sxy = (sx, sy): stretch figure from default size
    """

    ny, nx = np.shape(AA)
    figxy = (nx * sxy[0] / dpi, ny * sxy[1] / dpi)

    if ulimit:
        AA = np.clip(AA, ulimit[0], ulimit[1])

    nobox = (len(xypars) == 0)
    if nobox:
        xypars = (0, 0, 1, 1)
    x0, y0, dx, dy = xypars
    extent = (x0, x0 + dx*nx, y0, y0 + dy*ny)

    plt.figure(figname, figsize=figxy, dpi=dpi, tight_layout=True)
    plt.clf()
    plt.imshow(AA, cmap=cmap, aspect='auto', extent=extent)
    plt.title(caption)

    fig = plt.gca()
    if nobox:
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)

    plt_end(pause)


def graphB(B, figname='graphB', caption='', xpars=(0, 1), ulimit=(), sxy=(1, .3), line='#1f77b4', pause=0.):
    """
    xpars = (x0, dx); () no axes;
    sxy = (sx, sy): stretch figx by sx and figy by sy * figx
    """

    B = np.transpose(B)
    nx = len(B)
    sx, sy = sxy
    figxy = (nx * sx / dpi, nx * sx * sy / dpi)

    nobox = (len(xpars) == 0)
    if nobox:
        xpars = (0, 1)
    x0, dx = xpars
    xx = x0 + np.arange(nx) * dx

    plt.figure(figname, figsize=figxy, tight_layout=True)
    plt.clf()
    plt.plot(xx, B, line)

    plt.title(caption)
    plt.autoscale(enable=True, axis='x', tight=True)
    if ulimit:
        plt.ylim(ulimit)
    plt.grid(True)

    fig = plt.gca()
    if nobox:
        fig.xaxis.set_visible(False)

    plt_end(pause)


def plotAAB(AA, figname='plotAAB', capA='', capB='', cmap='gray', line='#1f77b4',
            xypars=(0, 0, 1, 1), ulimit=(), roi=(), sxy=(1, 1), aby=(3, 1), pause=0.):
    """
    xypars = (x0, y0, dx, dy); () no axes;
    sxy = (sx, sy): stretch figure from default size
    aby = (ay, by): relative y-size of AA and B
    """

    ''' plotAA '''
    ny, nx = np.shape(AA)
    sx, sy = sxy
    ay, by = aby
    figxy = (nx * sx / dpi, ny * sy * ((ay + by) / ay) / dpi)

    if ulimit:
        AA = np.clip(AA, ulimit[0], ulimit[1])

    nobox = (len(xypars) == 0)
    if nobox:
        xypars = (0, 0, 1, 1)
    x0, y0, dx, dy = xypars
    extent = (x0, x0 + dx*nx, y0, y0 + dy*ny)
    xx = x0 + np.arange(nx) * dx

    print(f'< plotAAB 1')
    plt.figure(figname, figsize=figxy, dpi=dpi, tight_layout=True)
    print(f'< plotAAB 2')
    plt.clf()
    plt.subplot2grid((ay + by, 1), (0, 0), rowspan=ay)
    plt.imshow(AA, cmap=cmap, aspect='auto', extent=extent)
    plt.title(capA)

    figa = plt.gca()
    if nobox:
        figa.axes.get_xaxis().set_visible(False)
        figa.axes.get_yaxis().set_visible(False)

    """ graphB """
    x1, x2, y1, y2 = extent
    if len(roi) == 0:
        roi = (nx//2, ny//2, 10, 10)
    rx0, ry0, rx, ry = roi
    iy = int(ny * (ry0 - y1) / (y2 - y1))
    B = AA[iy, :]
    put_cursor(roi, extent)

    plt.subplot2grid((ay + by, 1), (ay, 0), rowspan=by)
    plt.plot(xx, B, line)
    plt.title(capB)
    plt.autoscale(enable=True, axis='x', tight=True)
    if ulimit:
        plt.ylim(ulimit)
    plt.grid(True)

    figb = plt.gca()
    if nobox:
        figb.xaxis.set_visible(False)

    plt_end(pause)


def put_cursor(roi=(), extent=()):
    cline = 'yellow'
    croi = 'cyan'
    alpha = 0.75

    ax1, ax2, ay1, ay2 = extent
    if len(roi) == 0:
        roi = ((ax1 + ax2)/2, (ay1 + ay2)/2, (ax2 - ax1)/100, (ay2 - ay1)/100)
    rx0, ry0, rx, ry = roi

    rx1 = rx0 - rx/2
    rx2 = rx0 + rx/2
    ry1 = ry0 - ry/2
    ry2 = ry0 + ry/2

    plt.axhline(y=ry0, color=cline, alpha=alpha)
    # plt.axvline(x=rx0, color=cline, alpha=alpha)
    plt.plot([rx1, rx2, rx2, rx1, rx1], [ry1, ry1, ry2, ry2, ry1], color=croi, alpha=alpha)


def plt_end(pause):
    if pause:
        plt.pause(pause)
    else:
        plt.show()


def graph_many(graphs, figname='graph_many', col_row=(1, 1), xpars=(0, 1), sxy=(1, 1), line='#1f77b4', pause=0.):
    ncol, nrow = col_row
    sx, sy = sxy
    nx = len(graphs[0][0])
    figxy = (ncol * nx * sx / dpi, nrow * nx * sx * sy /dpi)

    nobox = (len(xpars) == 0)
    if nobox:
        xpars = (0, 1)
    x0, dx = xpars
    xx = x0 + np.arange(nx) * dx

    plt.figure(figname, figsize=figxy, tight_layout=True)
    plt.clf()

    for graph in graphs:
        uu, (irow, icol), caption, ylimit = graph

        index = (irow - 0) * ncol + (icol + 1)
        plt.subplot(nrow, ncol, index)
        plt.plot(xx, uu, line)
        plt.autoscale(enable=True, axis='x', tight=True)
        if ylimit:
            plt.ylim(ylimit)
        plt.title(caption)
        plt.grid(True)

        fig = plt.gca()
        if nobox:
            fig.xaxis.set_visible(False)

    plt_end(pause)


def mayaviAA(AA, figname='mayaviAAx', caption='', view=(70, 20), ulimit=(), sxy=(1, 1), cmap='jet'):
    import mayavi.mlab as ml

    el, az = view
    sx, sy = sxy
    figxy = (700 * sx, 500 * sy)

    if ulimit:
        AA = np.clip(AA, ulimit[0], ulimit[1])

    ml.figure(figname, size=figxy)
    ml.clf()
    ml.surf(AA, colormap=cmap, warp_scale='auto')
    ml.text(.02, .02, caption, width=len(caption)*.015)
    ml.view(elevation=el, azimuth=az)
    ml.show()


if __name__ == '__main__':
    pi2 = 2 * np.pi
    nx, ny = (1500, 1000)
    xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
    kk = pi2 / 500

    aa = 10 * np.sin(kk * xx) * np.sin(kk * yy) + 5 * (np.random.rand(1000, 1500) - .5)
    b = aa[125, :]

    graphs = []
    graphs += [(aa[0, :], (0, 0), 'aa[0, :], (0, 0)', ())]
    graphs += [(aa[1, :], (1, 0), 'aa[1, :], (1, 0)', ())]
    graphs += [(aa[2, :], (2, 0), 'aa[2, :], (2, 0)', ())]
    graphs += [(aa[3, :], (1, 1), 'aa[3, :], (1, 1)', ())]
    graphs += [(aa[4, :], (2, 2), 'aa[4, :], (2, 2)', ())]
    graphs += [(aa[5, :], (4, 2), 'aa[5, :], (4, 2)', ())]

    # plotAA(aa, figname='fig1', caption='AA', sxy=(.5, .5), cmap='jet', pause=1)
    # graphB(b, caption='B', xpars=(-100, .5), sxy=(.5, .5), pause=1)
    plotAAB(aa, ulimit=(), roi=(300, 375, 100, 100), sxy=(.5, .5), capA='AA', capB='B', pause=1)
    # graph_many(graphs, col_row=(3, 5), sxy=(.3, .3), xpars=(), pause=1)

    mayaviAA(aa, caption='asdfdsaadfadfasdfasdf123')


