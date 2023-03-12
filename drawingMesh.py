fig, ax = plt.subplots()
#ax.set_xscale('symlog', linthresh=0.1)
#ax.set_yscale('symlog', linthresh=0.1)

ax.set_adjustable("datalim")
ax.plot([0,0],[0,120],'o-r',dashes=[7,10],color='black')
ax.plot([0,120],[0,120],'o-r',dashes=[7,10],color='black')
ax.plot([120,0],[0,0],'o-r',dashes=[7,10],color='black')
ax.plot([0+displacements[0],0],[0+displacements[1],120],'o-r',color='r')
ax.plot([0+displacements[0],120],[0+displacements[1],120],'o-r',color='r')
ax.plot([120,0+displacements[0]],[0,0+displacements[1]],'o-r',color='r')

import matplotlib.ticker as ticker

#scale_x = 0.83333
#scale_y = 0.83333

#0.17280006912002762

#ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(120))
#ax.xaxis.set_major_formatter(ticks_x)

#ticks_y = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(120))
#ax.yaxis.set_major_formatter(ticks_y)

#plt.yticks([0,0.1,20,40,60,80,100,120],[0,0.1,20,40,60,80,100,120])	
#plt.xticks([0,0.1,20,40,60,80,100,120],[0,0.1,20,40,60,80,100,120])

#plt.yticks([0,0.1,10e-1,10e0,120],[0,0.1,10e-1,10e0,120])	
#plt.xticks([0,0.1,10e-1,10e0,120],[0,0.1,10e-1,10e0,120])
ax.set_aspect(1)


plt.show()