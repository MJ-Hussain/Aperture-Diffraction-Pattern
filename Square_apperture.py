
import numpy as np
import matplotlib.pyplot as plt


#------------------------------------------------------------#
#plot of Intensity
#------------------------------------------------------------#
def intensity_plot(x,y,z,S_id):
    X, Y = np.meshgrid(x, y)
    # cp = plt.contourf(X, Y, Intensity)
        
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, z, levels=100)
    ax.set_title('Differaction Pattern')
    ax.set_xlabel(r'sin$\theta$')
    ax.set_ylabel(r'sin$\Phi$')
    return
#------------------------------------------------------------#
#Goodness of fit test
#------------------------------------------------------------#

def chi_sq_test(Ob,Func):
    N=len(Ob)
    #n or constraints for our function
    n=2
    Sum_calc=0
    for i in range(N):
        Sum_calc+=(Ob[i]-Func[i])**2
    chi=Sum_calc/(N-n-1)
    return chi
#------------------------------------------------------------#
#Aperture dimension calculation function
#------------------------------------------------------------#
def dimension(dim_data,axis):
    length=(len(axis))
    def I_func(wi):
        
        I=np.zeros(length)
        for i in range(length):
            det=k*wi*0.5*axis[i]
            if det==0:
                I[i]=Io
            else:
                I[i]=Io*(np.sin(det)/det)**2
        return I        
    W=np.linspace(1,100,100)*1e-6
    
    chi=[]
    for wi in W:
        In=I_func(wi)        
        chi.append(chi_sq_test(dim_data,In))
    ind=chi.index(min(chi))
    req_dim=W[ind]
    I_fit=I_func(req_dim)
    return req_dim,min(chi),I_fit
#------------------------------------------------------------#


def ProcessData(file_name):
    
    start_read=0
    line_counter=0
    with open(file_name,'r') as file:
    
     for line in file:
        line_counter+=1
        line1=line.strip()
        if start_read==1:
            words=line1.split('=')
            if words[0]=='issid':
                issid=words[1]
            elif words[0]=='Wavelength (nm)':
                lemda=float(words[1])
            elif words[0]=='shape':
                shape=words[1]
            elif words[0]=='features':
                features=int(words[1])
            elif words[0]=='With Noise':
                Noise=words[1]
            elif words[0]=='Variable Wavelength':
                V_lemda=words[1]
            elif words[0]=='As Projected':
                projected=words[1]
            elif words[0]=='Horizontal, Vertical units':
                unit=words[1]
            
        if line1=='<MetaDataAtStart>':
            start_read=1
        elif line1=='</MetaDataAtStart>':
            break
    
    
    A=np.loadtxt(file_name, skiprows=line_counter+1, unpack=True)
    #Axis data
    Y_sin_phi=A[0,1:]
    X_sin_theta=A[1:,0]
    #Intensity data
    Intensity=A[1:,1:]
    

    

    
    #Slice of data from X and Y for function fit
    length=len(X_sin_theta)
    X_data=Intensity[:,int(length/2)]
    Y_data=Intensity[int(length/2),:]
    #Diagonal data
    diagonal_data=np.zeros(length)
    for i in range(length):
        diagonal_data[i]=Intensity[i,i]
    global Io,k
    #Max intensity refered to as Io
    Io=Intensity.max()
    #wave number k
    k=2*np.pi/lemda*1e9
    
        
    Wx,err_x,Ix=dimension(X_data,X_sin_theta)
    Ly,err_y,Iy=dimension(Y_data,Y_sin_phi)
    
    
    
    def xy_func_plot():
        
        plt.figure(2)
        plt.subplot(211)
        plt.plot(X_sin_theta, Ix, color='blue', label='Fitted Function')
        plt.plot(X_sin_theta, X_data, color='black', label=r'datapoints (Sin$\Phi$=0)')
        plt.xlabel(r'sin$\theta$')
        plt.ylabel(r'Intensity')
        plt.title(r'Fitted Function for points Sin$\Phi$=0')
        plt.legend()
        
    
        plt.subplot(212)
        plt.plot(Y_sin_phi, Iy, color='blue', label='Fitted Function')
        plt.plot(Y_sin_phi, Y_data, color='black', label=r'datapoints (Sin$\theta$=0)')
        plt.xlabel(r'sin$\Phi$')
        plt.ylabel(r'Intensity')
        plt.title(r'Fitted Function for points Sin$\theta$=0')
        plt.legend()
        
        plt.subplots_adjust(hspace=0.6)
        return
    
    area=Wx*Ly
    area_err=area*np.sqrt((err_x/Wx)**2+(err_y/Ly)**2)
    Io_per_A=Io/area
    error=Io*area_err/area**2
    
    
    #Call plot function
    intensity_plot(X_sin_theta,Y_sin_phi,Intensity,issid)
    xy_func_plot()
    
    
    results = {
        "shape": shape, # one of "square", "rectangle", "diamond", "circle" - must always be present.
        "dim_1": Wx, # a floating point number of the first dimension expressed in microns
        "dim_1_err": err_x, # The uncertainity in the above, also expressed in microns
        "dim_2": Ly, # For a rectangle, the second dimension, for other shapes, the same as dim_1
        "dim_2_err": err_y,  # The uncertainity in the above, also expressed in microns
        "I0/area": Io_per_A, # The fitted overall intensity value/area of the aperture.
        "I0/area_err": error, # The uncertainity in the above.
    }
    return results

if __name__=="__main__":
    from pprint import pprint
    filename='differaction_data.txt'
    test_results=ProcessData(filename)
    pprint(test_results)
    plt.show()
