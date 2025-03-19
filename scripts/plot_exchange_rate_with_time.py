# Author: Paul Guénon

import numpy as np
import argparse
import matplotlib.pyplot as plt

"""
This script reads a log file containing replica exchange data, calculates running averages of exchange rates, and plots the results over time.

Usage:
    python plot_exchange_rates_with_time.py -f <logfile> [-o <outputfile>] [-w <window>]

Arguments:
    -f, --file: Name of the log file (required)
    -o, --output: Name of the output PNG file (default: "exchange_rates_along_time.png")
    -w, --window: Window size in percentage for running average (default: 2)

Example:
    python plot_exchange_rates_with_time.py -f exchange_log.txt -o output.png -w 5
"""

def plot_progress_bar(current, total, bar_length=40):
    progress = current / total
    block = int(round(bar_length * progress))
    text = f"\rProgress: [{'#' * block + '-' * (bar_length - block)}] {progress * 100:.2f}%"
    print(text, end='')

def open_file (namefile) :
    file_opened=open(namefile,'r')
    lines_file=file_opened.readlines()
    data=[]
    count=0
    nreplica=0
    dt_tries=0
    times_restart=[]
    restart=False
    print("Loading file : "+namefile)
    for row in lines_file :
        plot_progress_bar(count, len(lines_file), bar_length=40)

        if "Repl ex" in row :
            data.append([x for x in row.split()])
            if nreplica==0 : 
                nreplica=int(row.split()[-1])+1

        if dt_tries==0 and "Replica exchange at step" in row :
            dt_tries=float(row.split()[-1])
        
        if "Restarting from checkpoint" in row : 
            restart=True
        if restart==True and "Replica exchange at step" in row : 
            times_restart.append(float(row.split()[-1]))
            restart=False
        
        count+=1
    plot_progress_bar(count, len(lines_file), bar_length=40)
    file_opened.close()
    return data, nreplica, dt_tries, times_restart

def running_average (y_data, window_size) :
    kernel = np.ones(window_size) / window_size
    y_avg = np.convolve(y_data, kernel, mode='valid')
    return y_avg

def main () :
    
    data,nreplica,dt_tries,times_restart=open_file(namefile)
    ntries=len(data)

    window_size=int(ntries*window_percentage/100)
    times_avg=np.arange((window_size-1)//2,(ntries+1)//2-(window_size)//2,1)*(dt_tries*2)/1000    

    color_choice=plt.cm.jet(np.linspace(0,1,nreplica-1))

    Exchange_array=np.zeros((nreplica-1,(ntries+1)//2))
    count_try=0
    print('\n')
    print("Processing data")
    for i in range (len(data)) :
        plot_progress_bar(i, len(data), bar_length=40)
        if len(data[i])>2 :
            if data[i][0]=='Repl' and data[i][1]=='ex' :
                for j in range(2,len(data[i])):
                    if data[i][j]=='x' :
                        Exchange_array[int(data[i][j-1]),int(count_try/2)]=1
                count_try+=1    
    plot_progress_bar(len(data), len(data), bar_length=40)
    print('\n')
    print("Running averages and plotting")

    ax=plt.subplots()[1]
    for c in range (nreplica-1) : 
        plot_progress_bar(c, nreplica-1, bar_length=40)
        ax.plot(times_avg,running_average(Exchange_array[c,:],window_size), color=color_choice[c],label=str(c)+'->'+str(c+1),lw=1)
    plot_progress_bar(nreplica-1, nreplica-1, bar_length=40)
    for time_restart in times_restart :
        ax.axvline(x=time_restart/1000, color='k', linestyle  = '--', lw=0.5)
    print('\n')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=4)
    ax.set_xlabel('Time in ns')
    ax.set_ylabel('Percentage of successfull exchanges')
    plt.title('Exchange rates from '+namefile,fontsize=8)
    plt.savefig(outputfile,dpi=300)
    plt.close()

if __name__ == "__main__":   
    parser = argparse.ArgumentParser(
        description="This script reads a log file containing replica exchange data, calculates running averages of exchange rates, and plots the results over time.",
        epilog="Example:\n    python plot_exchange_rates_with_time.py -f exchange_log.txt -o output.png -w 5",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    parser.add_argument('-f', '--file', required=True, help="Name of the log file (required)")
    parser.add_argument('-o', '--output', default="exchange_rates_along_time.png", help="Name of the output PNG file (default: 'exchange_rates_along_time.png')")
    parser.add_argument('-w', '--window', type=float, default=2, help="Window size in percentage of the trajectory for running average (default: 2)")

    args = parser.parse_args()
    namefile = args.file
    outputfile = args.output
    window_percentage = args.window
    main ()

