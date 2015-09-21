import numpy as num

def g(x,y,ma,mb,wx,wy,angle):
    
    rad=angle*num.pi/180.

    cos=num.cos(rad)
    coss=cos*cos
    sin=num.sin(rad)
    sins=sin*sin
    
    sin2r=num.sin(2*rad)
    wxs=wx*wx
    wys=wy*wy

    a=coss/(2*wxs)+sins/(2*wys)
    b=-sin2r/(4*wxs)+sin2r/(4*wys)
    c=sins/(2*wxs)+coss/(2*wys)


    return num.exp(-a*(x-ma)**2-2*b*(x-ma)*(y-mb)-c*(y-mb)**2)

def getKernel(wa,wb,ang):
    
    coords=[[0.5,0.5],[0.5,1.5],[0.5,2.5],[1.5,0.5],[1.5,1.5],[1.5,2.5],[2.5,0.5],[2.5,1.5],[2.5,2.5]]

    p=[]
    for vi in range(len(coords)):
        X=coords[vi][0]
        Y=coords[vi][1]

        p.append(g(X,Y,1.5,1.5,wa,wb,ang))

    norm=sum(p)

    kernel=num.array([[p[0],p[1],p[2]],[p[3],p[4],p[5]],[p[6],p[7],p[8]]])/norm

    return kernel

if __name__=='__main__':
    w=0.44
    print getKernel(w,w*2,0.)

