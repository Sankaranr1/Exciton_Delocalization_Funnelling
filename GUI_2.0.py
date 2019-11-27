# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 11:13:08 2019
@author: sankaran
"""
import time
import matplotlib
import cmath
from tkinter import *
from tkinter.ttk import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from scipy import *
from scipy.integrate import simps
from matplotlib import animation
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import math

"""Variables Initialization"""
n1=1
n2=1
b=1.0
tau=1.0
height=0.84
var_n3=0
var_n4=0

"""Tkinter window"""
window=Tk()
window.title("GUI")
window.geometry("1040x520")

"""Tk variables initialization"""
b=DoubleVar()
h=DoubleVar()
tau=DoubleVar()
well=IntVar()
vh=StringVar(value="1")
vb=StringVar(value="1")
vtau=StringVar(value="20")

tkvarn1 = IntVar(window)
tkvarn1.set(1)
tkvarn2 = IntVar(window)
tkvarn2.set(1)
tkvarn3 = IntVar(window)
tkvarn3.set(1)
tkvarn4 = IntVar(window)
tkvarn4.set(1)
choices1 = {0,1,2,3,4,5,6,7,8,9,10}
choices = {0,1,2,3,4,5,6,7,8,9,10}

"""Tk objects"""
label2=Label(window,text="Choose number of Wells:")
label2.place(relx=0.07,rely=0.03)

label1=Label(window,text="Enter the parameters") 
                              
labeln1=Label(window,text="n1:")                                                # well width n1
infn1=Label(window, text="Well 1 width (PbI6 units)")
popupMenun1 = OptionMenu(window, tkvarn1, *choices1)

labeln2=Label(window, text="n2:")                                               # well width n2
infn2=Label(window, text="Well 2 width (PbI6 units)")
popupMenun2 = OptionMenu(window, tkvarn2, *choices)

labelb=Label(window, text="b:")                                                 # barrier width b
infb=Label(window, text="Barrier width (nm)")
Eb=Entry(window,text="b:",textvariable=vb)

labelh=Label(window, text='V_0:')                                               # barrier height h
infh=Label(window, text="Height of Potential Barrier (V)")
Eh=Entry(window,text="V_0:",textvariable=vh)

labeltau=Label(window, text="\u03C4:")                                          # dephasing time constant tau
inftau=Label(window, text="Dephasing Time (fs)")
Etau=Entry(window,text="\u03C4:",textvariable=vtau)

labeln3=Label(window, text="n3:")                                               # well width n3
infn3=Label(window, text="Well 3 width (PbI6 units)")
popupMenun3 = OptionMenu(window, tkvarn3, *choices)
labeln4=Label(window, text="n4:")                                               # well width n4
infn4=Label(window, text="Well 4 width (PbI6 units)")
popupMenun4 = OptionMenu(window, tkvarn4, *choices)

figwell=Figure()                                                                #Figure of QW system for representational purpose
canvas = FigureCanvasTkAgg(figwell, master=window)
canvas.get_tk_widget().pack(side=LEFT)
#ax=figwell.add_subplot(111)
#ax.yaxis.set_major_locator(plt.NullLocator())
#ax.xaxis.set_major_formatter(plt.NullFormatter())

"""Start of function plotwell()


This function plots the correct representative image of the chosen QW system based on choice of user


"""
def plotwell():
    figwell.clf()
    ax=figwell.add_subplot(111)
    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    var=well.get()
    if (var==1):
        x=(0,1,1,1.63,1.63,2.63,2.63,3.26,3.26,4.26)
        y=(1,1,0,0,1,1,0,0,1,1)
    elif(var==2):
        x=(0,1,1,1.63,1.63,2.63,2.63,3.26,3.26,4.26,4.26,4.89,4.89,5.89)
        y=(1,1,0,0,1,1,0,0,1,1,0,0,1,1)
    else:
        x=(0,1,1,1.63,1.63,2.63,2.63,3.26,3.26,4.26,4.26,4.89,4.89,5.89,5.89,6.49,6.49,7.49)
        y=(1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1)
    ax.plot(x,y,color='red')
    canvas.draw()
    canvas.get_tk_widget().pack(side=LEFT)
#    toolbar = NavigationToolbar2Tk(canvas, window)
#    toolbar.update()
    MainWindow()

"""End of plotwell()"""

"""Start of function MainWindow

This function interfaces the user to input simulation parameters to the system

"""
def MainWindow():
    
    global var_n3,var_n4
    window.update()
    
    """Tk objects"""
    tkvarn1.set(1)        
    tkvarn2.set(1)        
    tkvarn3.set(1)        
    tkvarn4.set(1)        
    var=well.get()
    n_of_wells=var
    n_of_wells+=1
            
    label1.place(relx=0.58,rely=0.05)
    
    labelb.place(relx=0.53,rely=0.46)
    infb.place(relx=0.73,rely=0.46)
    Eb.place(relx=0.59,rely=0.46)
    
    labelh.place(relx=0.53,rely=0.53)
    infh.place(relx=0.73,rely=0.53)
    Eh.place(relx=0.59,rely=0.53)
    
    labeltau.place(relx=0.53,rely=0.60)
    inftau.place(relx=0.73,rely=0.60)
    Etau.place(relx=0.59,rely=0.60)

    labeln1.place(relx=0.53,rely=0.18)
    infn1.place(relx=0.73,rely=0.18)
    popupMenun1.place(relx=0.63,rely=0.18)
    labeln2.place(relx=0.53,rely=0.25)
    infn2.place(relx=0.73,rely=0.25)
    popupMenun2.place(relx=0.63,rely=0.25)

    if(n_of_wells==2):
        var_n3=1
        var_n4=1        
        
    if(n_of_wells==3):
        var_n3=0
        var_n4=1        
        labeln3.place(relx=0.53,rely=0.32)
        infn3.place(relx=0.73,rely=0.32)
        popupMenun3.place(relx=0.63,rely=0.32)
        
    if(n_of_wells==4):
        var_n3=0
        var_n4=0
        labeln3.place(relx=0.53,rely=0.32)
        infn3.place(relx=0.73,rely=0.32)
        popupMenun3.place(relx=0.63,rely=0.32)
        labeln4.place(relx=0.53,rely=0.39)
        infn4.place(relx=0.73,rely=0.39)
        popupMenun4.place(relx=0.63,rely=0.39)
        
    if(var_n4==1):
        popupMenun4.configure(state = DISABLED)
        tkvarn4.set(0)        
    else:
        popupMenun4.configure(state = NORMAL)        

    if(var_n3==1):
        popupMenun3.configure(state = DISABLED)
        tkvarn3.set(0)        
    else:
        popupMenun3.configure(state = NORMAL)        
    
    button1=Button(window, text="Run", command=display)
    button1.place(relx=0.7,rely=0.78)
    
   
"""End of MainWindow()"""

"""Start of function display()

This function runs the simulation and displays results

"""
def display():
    h=Eh.get()
    b=Eb.get()
    tau=Etau.get()
    n1=tkvarn1.get()
    n2=tkvarn2.get()
    n3=tkvarn3.get()
    n4=tkvarn4.get()
    
    print("values recorded")
    QWsimulation(n1,n2,n3,n4,b,tau,h)

"""End of display()"""    

"""Function:QWsimulation()

The calculating engine of the simulation

"""
def QWsimulation(n1,n2,n3,n4,b,tau,h):
    n_of_wells=well.get()
    n_of_wells+=1                                                               #Number of QWs
    
    """Variables"""
    n1=int(n1)
    n2=int(n2)
    n3=int(n3)
    n4=int(n4)
    tau=float(tau)
    height=float(h)
    b=float(b)
    b=b*1.0e-9 #in nm     

    a=zeros(n_of_wells)                                                         #width of wells
    a[0]=0.63e-9*n1 #in nm
    a[1]=0.63e-9*n2 #in nm

    if n_of_wells==3:
        a[2]=0.63e-9*n3 #in nm
    if n_of_wells==4:
        a[2]=0.63e-9*n3 #in nm        
        a[3]=0.63e-9*n4 #in nm

    Nx=zeros(4)
    if (n1!=0):Nx[0]=1
    if (n2!=0):Nx[1]=1
    if (n3!=0):Nx[2]=1
    if (n4!=0):Nx[3]=1

    N_wells=sum(Nx)
    leng=((sum(a))+(((1+N_wells)*b)))                                            #total length of the domain
            
    print("starting simulation")
    
    loading=Label(window, text="Evaluating result:")
    loading.place(relx=0.15,rely=0.86)
    progress=Progressbar(window,orient='horizontal',length=600,mode='determinate')
    progress.place(relx=0.25,rely=0.86)
    
    progress['value']=20                   
    window.update_idletasks()
        
    """Getting the correct ground state energy"""
    m=0.15                                                                      #effective mass of charge carrier
    m=m*(0.510998e6) / (9e16)                                                   #in eV s^2 m^-2
    h=6.582119e-16                                                              #hbar in eV.s (1.054e-34 is in J.s)
    #h=1.0546e-34
    C=2*m/(h**2)                                                                #in [m^-2 eV-1]
    #C= 5.236556965233972e+18
    N=1000                                                                      #Number of points in space
    VV=zeros(N)                                                                 #Defines potential at each point
    
    """interfaces"""   
    M=zeros(int(2*(N_wells)+2))                                                 #array M has coordinates of all the interfaces if wells and barriers
    ind=0

    for i in range(1,len(M)):
        if(i%2==1):
            M[i]=M[i-1]+b
        else:
            M[i]=M[i-1]+a[ind]
            ind+=1
            
    
    """function to calculate wavefunction""" 
    def FPW_demo(psi0,nu0,x,epsilon):        

        N=1000    
        xmin=x[0]
        xmax=x[N-1]
        dx=(xmax-xmin)/N
        dx2=dx*dx
        
        psi=zeros(N)                                                            #Wavefunction
        phi=zeros(N)
        
        psi[0]=psi0[0]                                                          #Initial condition
        psi[1]=psi0[1]
        
        def pot(x0):
            pe=1                                                                #defines the potential
            for i in range(1,len(M)):
                if (x0==0):
                    pe=nu0
                if (x0>M[i-1] and x0<M[i]):
                    if (mod(i,2)==0):
                        pe=0
                    else:
                        pe=nu0
            return pe
                   
        for i in range(1,N-1):
            x0=xmin+i*dx
            VV[i]=pot(x0)                                                       #ODE solver-Numerov's Method
            fact1=2.0*(1.0+((5.0*dx2/12.0)*(VV[i]-epsilon)))
            fact2=(1.0-((dx2/12.0)*(VV[i+1]-epsilon)))
            fact3=(1.0-((dx2/12.0)*(VV[i-1]-epsilon)))
            psi[i+1] = (fact1*psi[i] - fact3*psi[i-1])/fact2
            
        """normalising the wave function"""
        norm = simps(psi**2,x)
        psi *= 1/sqrt(norm)
        return (psi)
    
    """main documentation""" 
    progress['value']=50                   
    window.update_idletasks()

    Z=(n1+n2+n3+n4) + 2                                                         #upper bound for number of eigenfunctions
    nu0=height*C                                                                #height of potential well
    epsilon=0.1*C
    psi0=array([0.0,1.0])    
    x=linspace(0,leng,N)    
    window.update_idletasks()
    
    """finding eigenenergies: Shooting Method"""         
    energies=zeros(Z)                                                           #eigenenergies
    xmin=x[0]
    xmax=x[N-1]
    E1=0
    dE = 0.001*C
    tol=1e-5
    
    fun=zeros((Z,N))                                                            #Eigenstate wavefunctions
    fun2=zeros((Z,N))                                                           #Eigenstate probability density functions
        
    for j in range (Z):
        progress['value']=50+4*j
        window.update_idletasks()                    
        
        j=j+1
                                             
        """initial guess of energy bounds"""                
        psi=zeros(N)
        psi1 = FPW_demo(psi0,nu0,x,E1)      
        b1=psi1[N-1]
        count=0
        while True:
            E2=E1+dE
            psi2=FPW_demo(psi0,nu0,x,E2)
            b2=psi2[N-1]
            count=count+1
            if (count>5000):
                break
            if (b1*b2 < 0):
                break        
            else:
                E1=E2
                b1=b2
    
        psi1=FPW_demo(psi0,nu0,x,E1)
        psi2=FPW_demo(psi0,nu0,x,E2)
    
        psi1prob=zeros(N)
        psi2prob=zeros(N)
        
        for i in range(N):
            psi1prob[i]=psi1[i]*psi1[i]                                         #Probability Density function
            psi2prob[i]=psi2[i]*psi2[i]
        
        """bisection"""                              
        count=0
    
        def pe(x0):
            pet=1            
            for i in range(1,len(M)):
                if (x0==0):
                    pet=nu0
                
                if (x0>M[i-1] and x0<M[i]):
                    if (mod(i,2)==0):
                        pet=0
                    else:
                        pet=nu0
            return pet
        
        dx=(xmax-xmin)/N
        
        for i in range(N):
            x1 = xmin+i*dx
            VV[i]=pe(x1)
            
        x=linspace(0,leng,N)  
        run=0
        while True:
            run=run+1
            E0 = (E1+E2)/2
            psi1 = FPW_demo(psi0,nu0,x,E0)
            b0=psi1[N-1]
            err=abs(b0-0)
            if (count>100):
                break
            if (err < tol):
                break
            else:
                if (b0*b1 < 0):
                    E2=E0
                    b2=b0
                else:
                    E1=E0
                    b1=b0
        E0=(E1+E2)/2.0
    
        phi=FPW_demo(psi0,nu0,x,E0)                                             
        phiprob=zeros(N)
        for i in range(N):
            phiprob[i]=phi[i]*phi[i]
         
        fun[j-1]=phi                                                            #Eigenstate wavefunction
        fun2[j-1]=phiprob                                                       #Eigenstate probability density function
        
        energies[j-1]=E0
        E1=E0+0.001*C
      
    """
    Find and use only eigenfunctions which are bounded
    """      
    idx = where(energies < nu0 + 0.5)
    funs = fun[idx]
    Z = size(energies[idx])
    
    progress['value']=300                    
    window.update_idletasks()
    
    """Gaussian function"""    
    mu = b + a[0]/2                                                             #center of first well, for gaussian
    sigma=a[0]*0.25
    gauss=zeros(N)
    gauss=(1/(sigma*(sqrt(2*pi))))*exp(-((x-mu)**2/(2*(sigma**2))))
            
    """normalising the Gaussian"""
    norm=0.0
    norm = simps((gauss)**2, x)
    gauss*= 1/sqrt(norm)
    
    """inner product""" 
    c=zeros(Z)                                                                  #coefficients of the basis functions in the linear combination
    c2=zeros(Z)                                                                 #squared coefficients of the basis functions in the linear combination
    wf=zeros((N,Z))                                                             #Eigenstates (transposed)
    
    for i in range(Z):
        wf[:,i]=funs[i]                                                         #transpose of fun[]
        c[i]=dot(wf[:,i],gauss)                                                 
        c2[i]=c[i]*c[i]
        
    c_weights = zeros(Z)                                                        #finding weights-to normalize the coefficients
    c_sum = sqrt(sum(c2))
    
    for i in range(Z):
        c_weights[i] = c[i]/c_sum
        
    c=c_weights
    c2[:]=c[:]*c[:]    
    
    wfn=zeros(N)                                                                #Gaussian wavefunction as Linear combination of basis functions
    for i in range(Z):
        wfn+=(c[i]*funs[i])
            
    norm = simps(wfn**2, x)                                                     #normalising the gaussian wave function
    wfn *= 1/sqrt(norm)
    
    progress['value']=375
    window.update_idletasks()
    
    """
    Time evolution
    
    Version 2: damping the oscillatory part of probability density
    N = x steps
    Z = basis functions
    T = time
    """ 
    T = 200                                                                     #time window in fs
    t_steps=200                                                                 #Number of steps between 0 and T fs
    tim=linspace(0,T,t_steps)                
    E_n=zeros(Z)                                                                #energies of eigenstates, in eV
    E_n[:]=energies[0:Z]/C
    wfn0=zeros((N,Z),dtype=complex)                                             #complex wavefunction, component-wise
    wfn=zeros(N,dtype=complex)                                                  #total complex wavefunction
    wfntprob_c=zeros((N,Z,T),dtype=complex)                                     #Probability density function, component-wise, constant part 
    wfntprob_o=zeros((N,Z,T),dtype=complex)                                     #Probability density function, component-wise, oscillatory part 
    wfnt=zeros((N,T),dtype=complex)                                             #Total complex probability density function (oscillatory+constant)
    wfntprob=zeros((N,T))                                                       #FINAL Real probability density function
    
    for i in range(Z):
        
        wfn0[:,i]=c[i]*wf[:,i]                                                  #wavefunction of initial gaussian, t=0, in terms of each eigenstate: (c_k)*(phi_k)
        wfn[:]+=wfn0[:,i]
        
    """constant part of wavefunction"""
    for i in range(Z):
        wfntprob_c[:,i,0]+=(conj(wfn0[:,i])*wfn0[:,i])
        for g in range(1,T):     
            wfntprob_c[:,i,g]=wfntprob_c[:,i,0]
    
    
        
    """oscillatory part of wavefunction"""
    for g in range(T):
        for k in range(Z):
            for l in range(Z):
                if(k!=l):
                    wfntprob_o[:,k,g]+=(conj(wfn0[:,k])*wfn0[:,l])*exp(-1j*(E_n[l]-E_n[k])*(1/(h))*g*10**(-15))*exp(-g/tau)
    
    for g in range(T):
        for k in range(Z):
            wfnt[:,g]+=(wfntprob_c[:,k,g])
            wfnt[:,g]+=(wfntprob_o[:,k,g])
        
        """normalization"""    
        wfntprob = wfnt.real                                                    #real part of psi squared at each time step
        norm = simps(wfntprob[:,g], x)                                          #normalisation at each time step
        wfntprob[:,g]*= 1/(norm)

    progress['value']=450           
    window.update_idletasks()    
        
    """
    Calculating efficiency
    """    
    def efficiency(T, b, a, leng, psisq):
        global sumQR,sumST
#        global yQR,yUV,yST,yWX   
        sumQR=zeros((T,T))
        sumST=zeros((T,T))    
        lengQR=int((N/leng)*(b+a[0]))-int((N/leng)*b) 
        lengST=int((N/leng)*((2*b)+a[0]+a[1]))-int((N/leng)*((2*b)+a[0]))             
        yQR=zeros((lengQR,T))
        yST=zeros((lengST,T))    
        
        QR=linspace(b,b+a[0],lengQR)
        ST=linspace(2*b+a[0],2*b+a[0]+a[1],lengST)
        
        yQR[:,:]=psisq[int((N/leng)*b):int((N/leng)*(b+a[0])),:]
        yST[:,:]=psisq[int((N/leng)*((2*b)+a[0])):int((N/leng)*((2*b)+a[0]+a[1])),0:]
        
        sumQR[0,0]=simps(yQR[:,0],QR) 
        sumST[0,0]=simps(yST[:,0],ST)
        
        for u in range(1,T):
            for uu in range(u):
                sumQR[uu,u]=sumQR[uu,u-1] 
                sumST[uu,u]=sumST[uu,u-1] 
            sumQR[u,u]=simps(yQR[:,u],QR) 
            sumST[u,u]=simps(yST[:,u],ST) 
#        for u in range(T):            
#            sumQR[:,u]=sumQR[:,u]/(sumQR[0,u])
#            sumST[:,u]=sumST[:,u]/(sumQR[0,u])
            
        eff1=(sumST[0,0]/(sumQR[0,0]+sumST[0,0]))*100
        eff2=(sumST[(T-1),(T-1)]/(sumQR[(T-1),(T-1)]+sumST[(T-1),(T-1)]))*100
#        eff3=(sumST[0,0]/(sumQR[0,0]+sumST[0,0]))*100
#        eff4=(sumST[499,499]/(sumQR[499,499]+sumST[499,499]))*100
            
                
        if(n_of_wells==3):
            global sumUV
            sumUV=zeros((T,T))
            lengUV=int((N/leng)*(3*b+a[0]+a[1]+a[2]))-int((N/leng)*(3*b+a[0]+a[1]))
            yUV=zeros((lengUV,T))
            
            UV=linspace(3*b+a[0]+a[1],3*b+a[0]+a[1]+a[2],lengUV)
            yUV[:,:]=psisq[int((N/leng)*(3*b+a[0]+a[1])):int((N/leng)*(3*b+a[0]+a[1]+a[2])),:]            
            sumUV[0,0]=simps(yUV[:,0],UV) 
            for u in range(1,T):
                for uu in range(u):
                    sumUV[uu,u]=sumUV[uu,u-1] 
                sumUV[u,u]=simps(yUV[:,u],UV)
#            for u in range(T):            
#                sumUV[:,u]=sumUV[:,u]/(sumQR[0,u])    

        if(n_of_wells==4):
            sumUV=zeros((T,T))
            lengUV=int((N/leng)*(3*b+a[0]+a[1]+a[2]))-int((N/leng)*(3*b+a[0]+a[1]))
            yUV=zeros((lengUV,T))            

            UV=linspace(3*b+a[0]+a[1],3*b+a[0]+a[1]+a[2],lengUV)
            yUV[:,:]=psisq[int((N/leng)*(3*b+a[0]+a[1])):int((N/leng)*(3*b+a[0]+a[1]+a[2])),:]
            sumUV[0,0]=simps(yUV[:,0],UV) 
            for u in range(1,T):
                for uu in range(u):
                    sumUV[uu,u]=sumUV[uu,u-1] 
                sumUV[u,u]=simps(yUV[:,u],UV)       
                
            global sumWX
            sumWX=zeros((T,T))
            lengWX=int((N/leng)*(4*b+a[0]+a[1]+a[2]+a[3]))-int((N/leng)*(4*b+a[0]+a[1]+a[2]))
            yWX=zeros((lengWX,T))

            WX=linspace(4*b+a[0]+a[1]+a[2],4*b+a[0]+a[1]+a[2]+a[3],lengWX)
            yWX[:,:]=psisq[int((N/leng)*(4*b+a[0]+a[1]+a[2])):int((N/leng)*(4*b+a[0]+a[1]+a[2]+a[3])),:]
            sumWX[0,0]=simps(yWX[:,0],WX) 
            for u in range(1,T):
                for uu in range(u):
                    sumWX[uu,u]=sumWX[uu,u-1] 
                sumWX[u,u]=simps(yWX[:,u],WX) 
#            for u in range(T):            
#                sumUV[:,u]=sumUV[:,u]/(sumQR[0,u])
#                sumWX[:,u]=sumWX[:,u]/(sumQR[0,u])

        print('efficiency=',eff1,eff2)
        return eff1,eff2

    eff1,eff2 = efficiency(T, b, a, leng, wfntprob)
    
    """
    Time evolves psisq over time t, plots V, psisq vs x.
    """
    progress['value']=500
    window.update_idletasks()
    progress.destroy()
    loading.destroy()


    def time_evolution(x, psisq, V, t):
        fig=Figure()
        ax1=fig.add_subplot(211)
        ax2=fig.add_subplot(212)
        t_range=100
        line, = ax1.plot(x,psisq[:,0])
        line2, = ax2.plot(tim[0:t_range],sumQR[0:t_range,0],label="well 1")
        line3, = ax2.plot(tim[0:t_range],sumST[0:t_range,0],label="well 2")
        
        if(n_of_wells==3):
            line4, = ax2.plot(tim[0:t_range],sumUV[0:t_range,0],label="well 3")
        if(n_of_wells==4):
            line4, = ax2.plot(tim[0:t_range],sumUV[0:t_range,0],label="well 3")
            line5, = ax2.plot(tim[0:t_range],sumWX[0:t_range,0],label="well 4")

        
        ax2.legend(loc="upper right")
        ax2.set_xlabel("time(fs)")
        ax2.set_ylabel("probability")

        time_text = ax1.text(0.02, 0.95, '', fontsize=10, transform=ax1.transAxes)
        ax1.plot(x,(V/nu0)*max(psisq[:,0]))
        ax1.set_title('Transfer from n = ' + str(n1) + ' to n = ' + str(n2)+', b='+str(b*1e9)+' nm', fontsize=10)
#        ax1.text(0.02, 0.50, '$\eta = $' + str(round(eff2, 1)) + ' %', transform=ax1.transAxes)
        ax1.set_xlabel('x (nm)', fontsize=10)
        ax1.set_ylabel('$|\Psi(x)|^2$', fontsize=10)
        
        
        def init():
            line.set_data(x,psisq[:,0])
            line2.set_data(tim[0:t_range],sumQR[0:t_range,0])
            line3.set_data(tim[0:t_range],sumST[0:t_range,0])            
            time_text.set_text('')
            if(n_of_wells==3):
                line4.set_data(tim[0:t_range],sumUV[0:t_range,0])
                return line,line2,line3,line4,time_text
            elif(n_of_wells==4):
                line5.set_data(tim[0:t_range],sumWX[0:t_range,0])
                line4.set_data(tim[0:t_range],sumUV[0:t_range,0])
                return line,line2,line3,line4,line5,time_text            
            else:
                return line,line2,line3,time_text
        
        def animate1(i):
            line.set_ydata(psisq[:,i])                                          # update the data.
            line2.set_ydata(sumQR[0:t_range,i])  
            line3.set_ydata(sumST[0:t_range,i])  
            time_text.set_text('t = {0} fs'.format(i))
            if(n_of_wells==3):
                line4.set_ydata(sumUV[0:t_range,i])  
                return line,line2,line3,line4,time_text
            elif(n_of_wells==4):
                line5.set_ydata(sumWX[0:t_range,i])  
                line4.set_ydata(sumUV[0:t_range,i])  
                return line,line2,line3,line4,line5,time_text            
            else:
                return line,line2,line3,time_text
        
        timestr = time.strftime("%Y%m%d-%H%M")                                  #current date and time
        
        win2=Tk()
        win2.geometry("400x600")
        win2.title("Simulation window")
        l1=Label(win2,text="dynamic evolution:")
        l1.place(relx=0.1,rely=0.05)
        canvas = FigureCanvasTkAgg(fig, master=win2)
        canvas.draw()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        ani1 = animation.FuncAnimation(fig, animate1, frames=t_range, interval=50, blit=True, save_count=t)
        win2.mainloop()

    time_evolution(x, wfntprob, (VV), T)

"""End of QWsimulation()"""
wells2=Radiobutton(window,text="2 wells",variable=well,value=1,command=plotwell)    #Radiobutton to choose number of QWs
wells3=Radiobutton(window,text="3 wells",variable=well,value=2,command=plotwell)
wells4=Radiobutton(window,text="4 wells",variable=well,value=3,command=plotwell)

wells2.place(relx=0.07,rely=0.1)
wells3.place(relx=0.16,rely=0.1)
wells4.place(relx=0.25,rely=0.1)

window.mainloop()




