import numpy as np
import math as m
import matplotlib.pyplot as plt

def main():
    COES = np.loadtxt("data/COES.dat", delimiter=",", unpack=False)

    t=COES[:,0]

    COES = np.delete(COES, 0, 1)


    Dark=False

    plot_coes(t,COES,save_plot=True, title='System COEs',dark=Dark)

def plot_coes(ts,coes,dark=False,hours=False,days=False,show_plot=False,save_plot=False,title='COEs',figsize=(20,12),dpi=200):

    if dark:
        plt.style.use('dark_background')

    #create figure and axses
    fig,axs=plt.subplots(nrows=2,ncols=3,figsize=figsize)

    #figure titles
    fig.suptitle(title,fontsize=20)

    #x axis
    if hours: 
        ts=ts/3600
        xlabel='Time (hours)'
    elif days:
        ts=ts/3600/24
        xlabel='Time (days)'
    else:
        xlabel='Time (seconds)'


    #plot true anomaly
    axs[0,0].plot(ts,coes[:,3]*(360/(2*np.pi)))
    axs[0,0].set_title('True anomaly vs.Time', y=1.08)
    axs[0,0].grid(True)
    axs[0,0].set_ylabel('Angle (degrees)')

    #plot semi major axis
    axs[1,0].plot(ts,coes[:,0])
    axs[1,0].set_title('Semi Major Axis vs.Time', y=1.08)
    axs[1,0].grid(True)
    axs[1,0].set_ylabel('Semi Major Axis (km??)')
    axs[1,0].set_xlabel(xlabel)

    #plot eccentricity
    axs[0,1].plot(ts,coes[:,1])
    axs[0,1].set_title('Eccentricity vs.Time', y=1.08)
    axs[0,1].grid(True)

    #plot argument of periage
    axs[0,2].plot(ts,coes[:,4])
    axs[0,2].set_title('Argument of Periapse vs.Time',y=1.08)
    axs[0,2].grid(True)

    #plot inclination
    axs[1,1].plot(ts,coes[:,2]*(360/(2*np.pi)))
    axs[1,1].set_title('Inclination vs.Time', y=1.08)
    axs[1,1].grid(True)
    axs[1,1].set_ylabel('Angle (degrees)')
    axs[1,1].set_xlabel(xlabel)

    #plot RAAN
    axs[1,2].plot(ts,coes[:,5])
    axs[1,2].set_title('RAAN vs.Time',y=1.08)
    axs[1,2].grid(True)
    axs[1,2].set_xlabel(xlabel)

    plt.subplots_adjust(wspace=0.3)
    plt.subplots_adjust(hspace=0.4)

    if show_plot:
        plt.show()
    if save_plot:
        plt.savefig('results/'+title+'.png',dpi=dpi)


if __name__ == '__main__':
    main()