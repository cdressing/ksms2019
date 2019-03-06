import pylab as pl
def histograms(spf):
    bins = 20

    fig = pl.figure(figsize=(16,10.5))
    pl.subplots_adjust(hspace=0.6)

    pl.subplot(3,3,1)
    pl.hist(spf['Vmag'], bins=bins)
    pl.xlabel('Vmag')

    pl.subplot(3,3,2)
    pl.hist(spf['Ksig'], bins=bins*2)
    pl.xlabel('K/$\sigma_K$')
    pl.xlim(4, 13)

    pl.subplot(3,3,5)
    pl.hist(spf['Rs'], bins=bins)
    pl.xlabel('Rs')

    pl.subplot(3,3,4)
    pl.hist(spf['Teff'], bins=bins)
    pl.xlabel('Teff')

    pl.subplot(3,3,7)
    pl.hist(spf['exptime']/3600, bins=bins)
    pl.xlabel('exptime [h]')

    pl.subplot(3,3,8)
    pl.hist(spf['Kerr'], bins=bins)
    pl.xlabel('Kerr [m/s]')


    ax = pl.subplot(3,3,3)
    ax.axis('off')

    # pl.annotate('{} observations @ {}k or 45 min'.format(nobs, kexp), xy=(0.1, 0.95))
    # pl.annotate('Number of targets = {}'.format(num_targets), xy=(0.1, 0.8))
    # pl.annotate('Number of 10 hour nights = {:.0f}'.format(num_nights), xy=(0.1, 0.65))

    ax = pl.subplot(3,3,6)
    ax.axis('off')

    #pl.annotate('Filters:\n{}'.format(filters), xy=(0.1, 0.6))
    # pl.annotate('50 brightest\n+{} closest\n+all multis ({})\n+all USPs ({})'.format(len(close),
    #                                                                                  len(multi),
    #                                                                                  len(usp)),
    #             xy=(0.1, 0.0))

    #pl.savefig(outfile)


def scatter_plots(spf):
    fig = pl.figure(figsize=(16, 8.5))
    pl.subplots_adjust(hspace=0.5, wspace=0.4)

    pl.subplot(2, 3, 1)
    x = 'Vmag'
    y = 'Rp'
    pl.semilogy(spf[x], spf[y], 'ko')
    yt = pl.yticks()[0]
    pl.yticks(yt, yt)
    pl.xlabel('{}'.format(x))
    pl.ylabel('{}'.format(y))
    pl.ylim(1, 13)
    pl.xlim(7, 14)

    ax = pl.subplot(2, 3, 2)
    x = 'Perp'
    y = 'Rp'
    ax.loglog(spf[x], spf[y], 'ko')
    xt = pl.xticks()[0]
    pl.xticks(xt, xt)
    yt = pl.yticks()[0]
    pl.yticks(yt, yt)
    pl.xlabel('{}'.format(x))
    pl.ylabel('{}'.format(y))
    pl.xlim(0.3, 50)
    pl.ylim(1, 13)

    ax = pl.subplot(2, 3, 4)
    x = 'sinc'
    y = 'Rp'
    ax.loglog(spf[x], spf[y], 'ko')
    xt = pl.xticks()[0]
    pl.xticks(xt, xt)
    yt = pl.yticks()[0]
    pl.yticks(yt, yt)
    pl.xlabel('{}'.format(x))
    pl.ylabel('{}'.format(y))
    pl.xlim(2e4, 10)
    pl.ylim(1, 13)

    ax = pl.subplot(2, 3, 5)
    x = 'Vmag'
    y = 'Jmag'
    ax.plot(spf[x], spf[y], 'ko')
    pl.xlabel('{}'.format(x))
    pl.ylabel('{}'.format(y))
    # pl.xlim(-0.5, 0.5)
    # pl.ylim(1, 1.7)

    ax = pl.subplot(2, 3, 3)
    x = 'exptime'
    y = 'Rp'
    ax.semilogx(spf[x] / 3600, spf[y], 'ko')
    pl.xlabel('{} [h]'.format(x))
    pl.ylabel('{}'.format(y))
    pl.xticks([3, 10, 30, 100], [3, 10, 30, 100])
    # pl.xlim(-0.5, 0.5)
    # pl.ylim(1, 1.7)

    ax = pl.subplot(2, 3, 6)
    x = 'Kerr'
    y = 'Rp'
    ax.plot(spf[x], spf[y], 'ko')
    pl.xlabel('{} [m/s]'.format(x))
    pl.ylabel('{}'.format(y))
    # pl.xticks([3, 10, 30, 100], [3, 10, 30, 100])
    # pl.xlim(-0.5, 0.5)
    # pl.ylim(1, 1.7)
    #pl.savefig(outfile)


def hr_plot(spf):
    pl.semilogy(spf['Teff'], spf['lum'], 'ko')
    pl.xlabel('Teff [K]')
    pl.ylabel('L$_{\star}$ [L$_{\odot}$]')

    xl = pl.xlim()
    pl.xlim(xl[::-1])

    yt = pl.yticks()[0]
    pl.yticks(yt, yt)
    pl.ylim(0.001, 3)
