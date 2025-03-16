import jax.numpy as np
import matplotlib.pyplot as plt
from PIL import Image

import moving_least_squares_deformation as mlsd

def demo():
    image = np.array(Image.open("toy.jpg"))
    image = np.swapaxes(image,1,0)[:,::-1]

    width, height, _ = image.shape
    print(image.shape)
    vx,vy = np.meshgrid(np.arange(width),np.arange(height), indexing='ij')
    vx=vx.reshape(-1)
    vy=vy.reshape(-1)
    X = np.stack((vx,vy),axis=-1).astype(np.float64)

    P = np.array([
        [30, 160], [125, 160], [225, 160],
        [100, 80], [160, 80], [85, 20], [180, 22]
    ]).astype(np.float64)

    Q = np.array([
        [42, 140], [125, 160], [235, 185],
        [90, 80], [140, 80], [85, 20], [180, 20]
    ]).astype(np.float64)

    x = mlsd.mlsd_vec(X, P,Q)
    F = mlsd.mlsdDefGrad_vec(X, P,Q)



    plt.rcParams['lines.markersize'] = 0.1
    fig, ax = plt.subplots(1,5, sharey='all', figsize=(12,4))

    xrng = (-0.2*width,1.2*width)
    yrng = (-0.2*height,1.3*height)

    ax[0].set_ylim(yrng)
    ax[0].set_xlim(xrng)
    ax[0].set_aspect('equal')
    plt.sca(ax[0])
    plt.scatter(X[:,0], X[:,1], c=image.reshape(-1,3)/255.0)
    ax[0].set_title("Original Image")

    ax[1].set_ylim(yrng)
    ax[1].set_xlim(xrng)
    ax[1].set_aspect('equal')

    plt.sca(ax[1])
    plt.scatter(x[:,0], x[:,1], c=image.reshape((-1,3))/255.0)
    ax[1].set_title("Affine Deformation")

    ax[2].set_ylim(yrng)
    ax[2].set_xlim(xrng)
    ax[2].set_aspect('equal')

    plt.sca(ax[2])
    plt.scatter(vx, vy, c = x[:,0]-vx)
    ax[2].set_title("x-displacement Ref")

    ax[3].set_ylim(yrng)
    ax[3].set_xlim(xrng)
    ax[3].set_aspect('equal')

    plt.sca(ax[3])
    plt.scatter(x[:,0], x[:,1], c = x[:,0]-vx)
    ax[3].set_title("x-displacement Def")

    ax[4].set_ylim(yrng)
    ax[4].set_xlim(xrng)
    ax[4].set_aspect('equal')

    plt.sca(ax[4])
    plt.scatter(x[:,0], x[:,1], c = F[:,0,1])
    ax[4].set_title("F_{01} Def")

    plt.tight_layout(w_pad=0.1)
    plt.show()

demo()
