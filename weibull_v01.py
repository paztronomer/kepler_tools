import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

x = np.arange(1,100.)/50.

def weib(x,u,a,y):
    # http://www.itl.nist.gov/div898/handbook/eda/section3/eda3668.htm
    # x >= u, x is vector
    # y,a > 0
    # y: shape parameter
    # u: location parameter
    # a: scale parameter
    return (y/a) * ((x-u)/a)**(y-1) * np.exp(-((x-u)/a)**y)

def gauss(x,u,s):
    return 1/(s*np.sqrt(2*np.pi)) * np.exp((-(x-u)**2)/(2*s**2))

yy=[1.,2.,3.4,5.,10.] # shape
aa=[1.,2.5,5.,10.,20]
ss=[.25,.5,.75,1.,1.25]


fig=plt.figure(figsize=(8.5,8))
ax1=fig.add_subplot(221)
ax2=fig.add_subplot(222)
ax3=fig.add_subplot(212)


x=np.linspace(0.001,5,1000)
x2=np.linspace(-5,20,1000)
x3=np.linspace(0.001,40,1000)

for gamma in yy:
    w=weib(x=x,u=0,a=1.,y=gamma)
    ax1.plot(x,w,label='$\gamma$=%.1f'%(gamma),lw=2)

    #count, bins, ignored = plt.hist(np.random.weibull(5.,1000))
    #x = np.arange(1,100.)/50.
    #scale = count.max()/weib(x, 1., 5.).max()
    #plt.plot(x, weib(x, 1., 5.)*scale)

for alpha in aa:
    w2=weib(x=x3,u=0,a=alpha,y=5.)
    ax2.plot(x3,w2,label='$\eta$=%.1f'%(alpha),lw=2)

for sigma in ss:
    g=gauss(x=x2,u=5.,s=sigma)
    ax3.plot(x2,g,label='$\sigma$=%.1f'%(sigma),lw=2)


# legend
ax1.legend(frameon=False,loc='upper left')#(fancybox=False, framealpha=0.5)
ax2.legend(frameon=False)#(fancybox=True, framealpha=0.5)
ax3.legend(frameon=False)#(fancybox=True, framealpha=0.5)

# title
ax1.set_title('Weibull PDF, $\mu$=0, $\eta$=1',fontsize=13)
ax2.set_title('Weibull PDF, $\mu$=0, $\gamma$=5',fontsize=13)
ax3.set_title('Gaussian PDF, $\mu$=5',fontsize=13)

# labels
ax1.set_xlabel('x')
ax1.set_ylabel('Probability Density')
ax2.set_xlabel('x')
ax2.set_ylabel('Probability Density')
ax3.set_xlabel('x')
ax3.set_ylabel('Probability Density')

# limits
ax1.set_xlim([0.05,5.])
ax2.set_xlim([0.2,28])
ax3.set_xlim([2.,10.])

# scale
ax1.set_xscale('log')
ax2.set_xscale('log')
ax3.set_xscale('log')

# tics
ax3.set_xticks([2, 3, 4, 5, 6, 7, 8, 9, 10])
#ax3.set_xticklabels(["20", "200", "500"])
ax3.get_xaxis().set_major_formatter(ticker.ScalarFormatter())

# adjust
plt.subplots_adjust(left=0.07, top=.97, right=.985, bottom=.06, wspace=.20, hspace=.32)

plt.show()

fig.savefig('weibull.pdf', dpi=100, facecolor='w', edgecolor='w', \
            orientation='portrait', papertype='letter', format='pdf',\
            transparent=False, bbox_inches=None, pad_inches=0.1, \
            frameon=None)
