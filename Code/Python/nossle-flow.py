import jax.numpy as np
import numpy as onp
import matplotlib.pyplot as plt
from PIL import Image

import moving_least_squares_deformation as mlsd

def demo():
    plt.rcParams['lines.markersize'] = 0.5
    sz = 260
    sp = np.floor((sz/20)).astype(np.int32)
    print(sz)
    print(sp)
    print(sz/2)
    image = onp.zeros((sz,sz,3), dtype=np.int32)
    image[::sp,:,:] = 255
    image[:,::sp,:] = 255
    image = np.array(image)

    width, height, _ = image.shape
    print(image.shape)
    vx,vy = np.meshgrid(np.arange(width),np.arange(height), indexing='ij')
    vx=vx.reshape(-1)
    vy=vy.reshape(-1)
    X = np.stack((vx,vy),axis=-1).astype(np.float64)

    P = np.array([[j,0] for j in np.arange(0,sz+sp,sp)])
    P = np.concatenate((P,np.array([[j,sz] for j in np.arange(0,sz+sp,sp)])), axis=0)
    P = np.concatenate((P,np.array([[sz,j] for j in np.arange(0,sz+sp,sp)])), axis=0)
    P = np.concatenate((P,np.array([[0,j] for j in np.arange(0,sz+sp,sp)])), axis=0).astype(np.float64)
    #Q = np.array([[j,0+0.625*sz*(j/(sz))] for j in np.arange(0,sz+sp,sp)])
    #Q = np.concatenate((Q,np.array([[j,sz-0.25*sz*(j/sz)] for j in np.arange(0,sz+sp,sp)])), axis=0)
    #Q = np.concatenate((Q,np.array([[sz,0.625*sz+0.125*sz*(j/sz)] for j in np.arange(0,sz+sp,sp)])), axis=0).astype(np.float64)
    #
    Q = np.array([[j,0+0.375*sz*(j/(sz/2))] for j in np.arange(0,sz/2+sp,sp)])
    Q = np.concatenate((Q,np.array([[j,0.375*sz] for j in np.arange(sz/2+sp,sz+sp,sp)])), axis=0)
    Q = np.concatenate((Q,np.array([[j,sz-0.25*sz*(j/sz)] for j in np.arange(0,sz+sp,sp)])), axis=0)
    Q = np.concatenate((Q,np.array([[sz+5*sp*(0.25-(j/sz-0.5)**2),0.375*sz+0.375*sz*(j/sz)] for j in np.arange(0,sz+sp,sp)])), axis=0)
    Q = np.concatenate((Q,np.array([[0+5*sp*(0.25-(j/sz-0.5)**2),j] for j in np.arange(0,sz+sp,sp)])), axis=0).astype(np.float64)

    x = mlsd.mlsd_vec(X, P,Q)
    F = mlsd.mlsdDefGrad_vec(X, P,Q)
    fig, ax = plt.subplots(1,2, sharey='all', figsize=(10,5))

    xrng = (-0.2*width,1.2*width)
    yrng = (-0.2*height,1.3*height)

    ax[0].set_ylim(yrng)
    ax[0].set_xlim(xrng)
    ax[0].set_aspect('equal')
    plt.sca(ax[0])
    plt.scatter(X[:,0], X[:,1], c=image.reshape(-1,3)/255.0)
    plt.colorbar()
    ax[0].set_title("original image")

    ax[1].set_ylim(yrng)
    ax[1].set_xlim(xrng)
    ax[1].set_aspect('equal')

    plt.sca(ax[1])
    #plt.scatter(x[:,0], x[:,1], c=image.reshape((-1,3))/255.0)
    #plt.scatter(x[:,0], x[:,1], c=F[:,1,0])
    plt.scatter(x[:,0], x[:,1], c=(x[:,0]-X[:,0]))
    plt.colorbar()
    #plt.scatter(Q[:,0], Q[:,1], c=(1.0,0.0,0.0))
    ax[1].set_title("horizontal displacement")

    #ax[2].set_ylim(yrng)
    #ax[2].set_xlim(xrng)
    #ax[2].set_aspect('equal')
    #
    #plt.sca(ax[2])
    #plt.scatter(vx, vy, c = x[:,0]-vx)
    #ax[2].set_title("x-displacement Ref")
    #
    #ax[3].set_ylim(yrng)
    #ax[3].set_xlim(xrng)
    #ax[3].set_aspect('equal')
    #
    #plt.sca(ax[3])
    #plt.scatter(x[:,0], x[:,1], c = x[:,0]-vx)
    #ax[3].set_title("x-displacement Def")
    #
    #ax[4].set_ylim(yrng)
    #ax[4].set_xlim(xrng)
    #ax[4].set_aspect('equal')
    #
    #plt.sca(ax[4])
    #plt.scatter(x[:,0], x[:,1], c = F[:,0,1])
    #ax[4].set_title("F_{01} Def")

    plt.tight_layout(w_pad=0.1)
    plt.show()

demo()
