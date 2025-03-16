import jax.numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.animation as ani
from PIL import Image

import moving_least_squares_deformation as mlsd
from jax import vmap

def load_image(fn):
    image = np.array(Image.open("toy.jpg"))
    image = np.swapaxes(image,1,0)[:,::-1]

    width, height, _ = image.shape
    print(image.shape)
    vx,vy = np.meshgrid(np.arange(width),np.arange(height), indexing='ij')
    vx=vx.reshape(-1)
    vy=vy.reshape(-1)
    X = np.stack((vx,vy),axis=-1).astype(np.float64)

    return X, image, width, height

def generate_control_point_snapshots():
    # Define control points
    #
    # Add 0.5's to avoid NaNs in defgrad, etc
    P = np.array([
        [50, 180],
        [145, 180],
        [245, 180],
        [120, 100],
        [180, 100],
        [105, 40],
        [200, 42]
    ]).astype(np.float64) + 0.5

    Q = np.array([
        [32, 210],
        [145, 180],
        [225, 145],
        [70, 100],
        [140, 100],
        [105, 40],
        [200, 42]
    ]).astype(np.float64) + 0.5

    QQ = np.array([
        [72, 150],
        [145, 180],
        [265, 215],
        [160, 100],
        [230, 100],
        [105, 40],
        [200, 42]
    ]).astype(np.float64) + 0.5

    # Define times for which motion control points specified: (nt)
    tt = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])

    def rigid(S, theta, translation):
        # rotate about center and translate
        center = np.array([167.5, 184.5])
        return ((S - center) @ np.array([[np.cos(theta), np.sin(theta)],[-np.sin(theta), np.cos(theta)]]) + center) + translation
    # Motion control points: (nt, ncp, dim)
    SS = np.stack((rigid( 1.00*P , 0.00, np.array([  0,   0])),
                   rigid( 1.06*Q , 0.40, np.array([150,  75])),
                   rigid( 1.12*P , 0.80, np.array([300, 150])),
                   rigid( 1.18*QQ, 1.20, np.array([450,  75])),
                   rigid( 1.25*P , 1.57, np.array([600,   0])),
                   rigid( 1.18*QQ, 1.20, np.array([450,- 75])),
                   rigid( 1.12*P , 0.80, np.array([300,-150])),
                   rigid( 1.06*Q , 0.40, np.array([150,- 75])),
                   rigid( 1.00*P , 0.00, np.array([000,   0]))
                   )
                  )

    deform_coeffs = mlsd.motion_2d_splines(tt, SS)

    return tt, SS, deform_coeffs


def demo(fn):
    X, image, width, height = load_image(fn)

    tt, SS, deform_coeffs = generate_control_point_snapshots()

    # plot ranges to be used with ref config
    xrng1 = (-0.5*width,1.3*width)
    #yrng1 = (-0.3*height,1.3*height)
    yrng1 = (-0.8*height,1.8*height)

    # plot ranges to be used for full motion
    xrng2 = (-0.5*width,3.0*width)
    yrng2 = (-0.8*height,1.8*height)


    frames = 40


    def eulerian_view():
        #
        # Eulerian view of motion and path/streaklines
        #

        plt.rcParams['lines.markersize'] = 0.1
        fig = plt.figure()
        ims = []

        Xpath_pt = np.array([300.0,100.0])
        xstreak_pt = np.array([300.0,100.0])

        global Xstreak
        Xstreak = np.array([xstreak_pt])

        print("-- generating eulerian view --")

        plt.gca().set_xlim(xrng2)
        plt.gca().set_ylim(yrng2)
        plt.gca().set_aspect('equal')
        # plot pathline and streakline
        xp, = plt.plot([], [], 'r-')
        xs, = plt.plot([], [], 'b-')
        # plot configuration
        scat = plt.scatter(X[:,0], X[:,1], c=image.reshape(-1,3)/255.0)

        def animate(frame):
            global Xstreak

            print("Frame {:3}".format(frame))
            t = frame*(tt[-1]/frames)

            # compute motion
            x = mlsd.mlsd_2d_motion_vec(X, t, tt, deform_coeffs)

            # compute location of pathline particle
            xpath_pt = mlsd.mlsd_2d_motion_onepoint(Xpath_pt, t, tt, deform_coeffs)
            xpx,xpy = xp.get_data()
            xpx = np.append(xpx, xpath_pt[0])
            xpy = np.append(xpy, xpath_pt[1])
            xp.set_data(xpx,xpy)

            # compute reference particle at streakline position
            Xstreak_pt = mlsd.mlsd_2d_inv_motion_onepoint(xstreak_pt, t, tt, deform_coeffs)
            Xstreak = np.append(Xstreak, np.array([Xstreak_pt]), axis=0)
            # convect streakline to current configuration
            xstreak = np.swapaxes(mlsd.mlsd_2d_motion_vec(Xstreak, t, tt, deform_coeffs), 0, 1)
            xs.set_data(xstreak)

            scat.set_offsets(x)

        mov = ani.FuncAnimation(fig, animate, frames=frames, repeat=True)
        #plt.show()
        mov.save(filename="motion-eularian-view.mp4", writer ="ffmpeg")


    def eulerian_vel_quiver_view(pltimg=False):
        #
        # Eulerian view of velocity quiver plot
        #

        plt.rcParams['lines.markersize'] = 0.1
        fig = plt.figure()
        ims = []

        xx = np.arange(xrng2[0],xrng2[1],(xrng2[1]-xrng2[0])/15.0)
        yy = np.arange(yrng2[0],yrng2[1],(yrng2[1]-yrng2[0])/15.0)

        xxg,yyg = np.meshgrid(xx,yy)
        xxg=xxg.reshape(-1)
        yyg=yyg.reshape(-1)
        x = np.stack((xxg,yyg), axis=-1)

        print("-- generating eulerian vel quiver view --")
        # compute motion
        v = mlsd.mlsd_2d_motion_vel_vec(x, 0.0, tt, deform_coeffs)

        plt.gca().set_xlim(xrng2)
        plt.gca().set_ylim(yrng2)
        plt.gca().set_aspect('equal')
        if pltimg:
            y = mlsd.mlsd_2d_motion_vec(X, 0.0, tt, deform_coeffs)
            scat = plt.scatter(y[:,0], y[:,1], c=image.reshape(-1,3)/255.0)
        quiv =  plt.quiver(xxg,yyg, v[:,0], v[:,1], scale=15000, width=0.002)

        def animate(frame):
            t = frame*(tt[-1]/frames)
            if pltimg:
                y = mlsd.mlsd_2d_motion_vec(X, t, tt, deform_coeffs)
                scat.set_offsets(y)
            v = mlsd.mlsd_2d_motion_vel_vec(x, t, tt, deform_coeffs)
            quiv.set_UVC(v[:,0], v[:,1])

        mov = ani.FuncAnimation(fig, animate, frames=frames, repeat=True)
        #plt.show()
        if pltimg:
            mov.save(filename="motion-eularian-vel-quiver-with-image-view.mp4", writer ="ffmpeg")
        else:
            mov.save(filename="motion-eularian-vel-quiver-view.mp4", writer ="ffmpeg")


    def lagrangian_and_eularian_view():
        #
        # Lagrangian and Eularian view with path/streaklines
        #

        plt.rcParams['lines.markersize'] = 0.1
        fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw={'width_ratios': [1, (xrng2[1]-xrng2[0])/(xrng1[1]-xrng1[0])]})
        plt.tight_layout()

        Xpath_pt = np.array([300,100]).astype(np.float64)
        xstreak_pt = np.array([300,100]).astype(np.float64)

        ax1.set_xlim(xrng1)
        ax1.set_ylim(yrng1)
        ax1.set_aspect('equal')
        scat1 = ax1.scatter(X[:,0], X[:,1], c=image.reshape(-1,3)/255.0)
        xs, = ax1.plot([], [], 'b-')

        ax2.set_xlim(xrng2)
        ax2.set_ylim(yrng2)
        ax2.set_aspect('equal')
        scat2 = ax2.scatter(X[:,0], X[:,1], c=image.reshape(-1,3)/255.0)
        xp, = ax2.plot([], [], 'r-')

        print("-- generating lagrangian and eulerian view --")

        def animate(frame):
            print("Frame {:3}".format(frame))
            t = frame*(tt[-1]/frames)

            # compute motion
            x = mlsd.mlsd_2d_motion_vec(X, t, tt, deform_coeffs)

            # compute location of pathline particle
            xpath_pt = mlsd.mlsd_2d_motion_onepoint(Xpath_pt, t, tt, deform_coeffs)
            xpx,xpy = xp.get_data()
            xpx = np.append(xpx, xpath_pt[0])
            xpy = np.append(xpy, xpath_pt[1])
            xp.set_data(xpx,xpy)

            # compute reference particle at streakline position
            Xstreak_pt = mlsd.mlsd_2d_inv_motion_onepoint(xstreak_pt, t, tt, deform_coeffs)
            xsx,xsy = xs.get_data()
            xsx = np.append(xsx, Xstreak_pt[0])
            xsy = np.append(xsy, Xstreak_pt[1])
            xs.set_data(xsx,xsy)

            # plot cur configuration
            scat2.set_offsets(x)

        mov = ani.FuncAnimation(fig, animate, frames=frames, repeat=True)
        #plt.show()
        mov.save(filename="motion-lagrangian-and-eulerian-view.mp4", writer ="ffmpeg")

    def velocity_field_view():
        #
        # velocity field components on Lagrangian and Eularian view
        #

        plt.rcParams['lines.markersize'] = 0.1
        fig = plt.figure()
        ax1_1 = plt.subplot(2,2,1)
        ax2_1 = plt.subplot(2,2,3)
        ax1_2 = plt.subplot(2,2,2)
        ax2_2 = plt.subplot(2,2,4)
        plt.tight_layout()

        vmin = -225
        vmax = 225
        cmap=mpl.colormaps['coolwarm']
        cmap.set_over('yellow')
        cmap.set_under('yellow')
        cmap.set_bad('cyan')
        nm=colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

        fig.colorbar(mpl.cm.ScalarMappable(norm=nm, cmap=cmap), cax=None, ax=ax2_1)
        fig.colorbar(mpl.cm.ScalarMappable(norm=nm, cmap=cmap), cax=None, ax=ax1_1).remove()

        # compute motion
        x = mlsd.mlsd_2d_motion_vec(X, 0.0, tt, deform_coeffs)
        v = mlsd.mlsd_2d_motion_vel_vec(x, 0.0, tt, deform_coeffs)

        # plot images
        ax1_1.set_xlim(xrng1)
        ax1_1.set_ylim(yrng1)
        ax1_1.set_aspect('equal')
        # plot ref configuration
        scat1 = ax1_1.scatter(X[:,0], X[:,1], c=image.reshape(-1,3)/255.0)

        ax1_2.set_xlim(xrng2)
        ax1_2.set_ylim(yrng2)
        ax1_2.set_aspect('equal')
        # plot cur configuration
        scat2 = ax1_2.scatter(x[:,0], x[:,1], c=image.reshape(-1,3)/255.0)

        # plot x-component of vel
        ax2_1.set_title("v_0")
        ax2_1.set_xlim(xrng1)
        ax2_1.set_ylim(yrng1)
        ax2_1.set_aspect('equal')
        # plot ref configuration
        ax2_1.scatter(X[:,0], X[:,1], c=v[:,0], cmap=cmap, norm=nm)

        ax2_2.set_title("v_0")
        ax2_2.set_xlim(xrng2)
        ax2_2.set_ylim(yrng2)
        ax2_2.set_aspect('equal')
        # plot cur configuration
        ax2_2.scatter(x[:,0], x[:,1], c=v[:,0], cmap=cmap, norm=nm)

        print("-- generating velocity field view --")
        print("{:5}, {:>15}, {:>15}".format("Frame", "v_0 min", "v_0 max"))

        def animate(frame):
            t = frame*(tt[-1]/frames)

            print("{:5},".format(frame), end=" ")
            # compute motion
            x = mlsd.mlsd_2d_motion_vec(X, t, tt, deform_coeffs)
            v = mlsd.mlsd_2d_motion_vel_vec(x, t, tt, deform_coeffs)
            print("{:15.10f}, {:15.10f}".format(np.min(v[:,0]), np.max(v[:,0])))

            scat2.set_offsets(x)

            for collection in ax2_1.collections: collection.remove()
            ax2_1.scatter(X[:,0], X[:,1], c=v[:,0], cmap=cmap, norm=nm)

            for collection in ax2_2.collections: collection.remove()
            ax2_2.scatter(x[:,0], x[:,1], c=v[:,0], cmap=cmap, norm=nm)

        mov = ani.FuncAnimation(fig, animate, frames=frames, repeat=True)
        #plt.show()
        mov.save(filename="motion-vel-field-view.mp4", writer ="ffmpeg")

    def lagrangian_strain_field_view():
        #
        # lagrangian strain field component on Lagrangian and Eularian view
        #
        def LagStrain(F):
            return 0.5 * (F.T @ F - np.eye(2))
        LagStrain_vec = vmap(LagStrain, in_axes=0, out_axes=0)

        plt.rcParams['lines.markersize'] = 0.1
        fig = plt.figure()
        ax1_1 = plt.subplot(2,2,1)
        ax2_1 = plt.subplot(2,2,3)
        ax1_2 = plt.subplot(2,2,2)
        ax2_2 = plt.subplot(2,2,4)
        plt.tight_layout()

        vmin = -0.35
        vmax = 1
        cmap=mpl.colormaps['coolwarm']
        cmap.set_over('yellow')
        cmap.set_under('yellow')
        cmap.set_bad('cyan')
        nm=colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

        fig.colorbar(mpl.cm.ScalarMappable(norm=nm, cmap=cmap), cax=None, ax=ax2_1)
        fig.colorbar(mpl.cm.ScalarMappable(norm=nm, cmap=cmap), cax=None, ax=ax1_1).remove()

        # compute motion
        x = mlsd.mlsd_2d_motion_vec(X, 0.0, tt, deform_coeffs)

        F = mlsd.mlsd_2d_motion_defgrad_vec(X, 0.0, tt, deform_coeffs)
        E = LagStrain_vec(F)

        # plot images
        ax1_1.set_xlim(xrng1)
        ax1_1.set_ylim(yrng1)
        ax1_1.set_aspect('equal')
        # plot ref configuration
        scat1 = ax1_1.scatter(X[:,0], X[:,1], c=image.reshape(-1,3)/255.0)

        ax1_2.set_xlim(xrng2)
        ax1_2.set_ylim(yrng2)
        ax1_2.set_aspect('equal')
        # plot cur configuration
        scat2 = ax1_2.scatter(x[:,0], x[:,1], c=image.reshape(-1,3)/255.0)

        # plot 00-component of E
        ax2_1.set_title("E_00")
        ax2_1.set_xlim(xrng1)
        ax2_1.set_ylim(yrng1)
        ax2_1.set_aspect('equal')
        # plot ref configuration
        ax2_1.scatter(X[:,0], X[:,1], c=E[:,0,0], cmap=cmap, norm=nm)

        ax2_2.set_title("E_00")
        ax2_2.set_xlim(xrng2)
        ax2_2.set_ylim(yrng2)
        ax2_2.set_aspect('equal')
        # plot cur configuration
        ax2_2.scatter(x[:,0], x[:,1], c=E[:,0,0], cmap=cmap, norm=nm)

        print("-- generating lagrangian strain field view --")
        print("{:5}, {:>15}, {:>15}".format("Frame", "E_00 min", "E_00 max"))

        def animate(frame):
            t = frame*(tt[-1]/frames)

            print("{:5},".format(frame), end=" ")
            # compute motion
            x = mlsd.mlsd_2d_motion_vec(X, t, tt, deform_coeffs)

            F = mlsd.mlsd_2d_motion_defgrad_vec(X, t, tt, deform_coeffs)
            E = LagStrain_vec(F)
            print("{:15.10f}, {:15.10f}".format(np.min(E[:,0,0]), np.max(E[:,0,0])))

            scat2.set_offsets(x)

            for collection in ax2_1.collections: collection.remove()
            ax2_1.scatter(X[:,0], X[:,1], c=E[:,0,0], cmap=cmap, norm=nm)

            for collection in ax2_2.collections: collection.remove()
            ax2_2.scatter(x[:,0], x[:,1], c=E[:,0,0], cmap=cmap, norm=nm)

        mov = ani.FuncAnimation(fig, animate, frames=frames, repeat=True)
        # plt.show()
        mov.save(filename="motion-lagrangian-strain-field-view.mp4", writer ="ffmpeg")

    def strain_field_view():
        #
        # strain field component on Lagrangian and Eularian view
        #
        def LagStrain(F):
            return 0.5 * (F.T @ F - np.eye(2))
        LagStrain_vec = vmap(LagStrain, in_axes=0, out_axes=0)
        def EulStrain(F):
            Finv = np.linalg.inv(F)
            return 0.5 * (np.eye(2) - Finv.T @ Finv)
        EulStrain_vec = vmap(EulStrain, in_axes=0, out_axes=0)

        plt.rcParams['lines.markersize'] = 0.1
        fig = plt.figure()
        ax1_1 = plt.subplot(2,2,1)
        ax2_1 = plt.subplot(2,2,3)
        ax1_2 = plt.subplot(2,2,2)
        ax2_2 = plt.subplot(2,2,4)
        plt.tight_layout()

        vmin = -0.35
        vmax = 1
        cmap=mpl.colormaps['coolwarm']
        cmap.set_over('yellow')
        cmap.set_under('yellow')
        cmap.set_bad('cyan')
        nm=colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

        fig.colorbar(mpl.cm.ScalarMappable(norm=nm, cmap=cmap), cax=None, ax=ax2_1)
        fig.colorbar(mpl.cm.ScalarMappable(norm=nm, cmap=cmap), cax=None, ax=ax1_1).remove()

        # compute motion
        x = mlsd.mlsd_2d_motion_vec(X, 0.0, tt, deform_coeffs)

        F = mlsd.mlsd_2d_motion_defgrad_vec(X, 0.0, tt, deform_coeffs)
        E = LagStrain_vec(F)
        e = EulStrain_vec(F)

        # plot images
        ax1_1.set_xlim(xrng1)
        ax1_1.set_ylim(yrng1)
        ax1_1.set_aspect('equal')
        # plot ref configuration
        scat1 = ax1_1.scatter(X[:,0], X[:,1], c=image.reshape(-1,3)/255.0)

        ax1_2.set_xlim(xrng2)
        ax1_2.set_ylim(yrng2)
        ax1_2.set_aspect('equal')
        # plot cur configuration
        scat2 = ax1_2.scatter(x[:,0], x[:,1], c=image.reshape(-1,3)/255.0)

        # plot 00-component of E
        ax2_1.set_title("E_00")
        ax2_1.set_xlim(xrng1)
        ax2_1.set_ylim(yrng1)
        ax2_1.set_aspect('equal')
        # plot ref configuration
        ax2_1.scatter(X[:,0], X[:,1], c=E[:,0,0], cmap=cmap, norm=nm)

        ax2_2.set_title("e_00")
        ax2_2.set_xlim(xrng2)
        ax2_2.set_ylim(yrng2)
        ax2_2.set_aspect('equal')
        # plot cur configuration
        ax2_2.scatter(x[:,0], x[:,1], c=e[:,0,0], cmap=cmap, norm=nm)

        print("-- generating lagrangian strain field view --")
        print("{:5}, {:>15}, {:>15}".format("Frame", "E_00 min", "E_00 max"))

        def animate(frame):
            t = frame*(tt[-1]/frames)

            print("{:5},".format(frame), end=" ")
            # compute motion
            x = mlsd.mlsd_2d_motion_vec(X, t, tt, deform_coeffs)

            F = mlsd.mlsd_2d_motion_defgrad_vec(X, t, tt, deform_coeffs)
            E = LagStrain_vec(F)
            e = EulStrain_vec(F)
            print("{:15.10f}, {:15.10f}".format(np.min(E[:,0,0]), np.max(E[:,0,0])))

            scat2.set_offsets(x)

            for collection in ax2_1.collections: collection.remove()
            ax2_1.scatter(X[:,0], X[:,1], c=E[:,0,0], cmap=cmap, norm=nm)

            for collection in ax2_2.collections: collection.remove()
            ax2_2.scatter(x[:,0], x[:,1], c=e[:,0,0], cmap=cmap, norm=nm)

        mov = ani.FuncAnimation(fig, animate, frames=frames, repeat=True)
        # plt.show()
        mov.save(filename="motion-strain-field-view.mp4", writer ="ffmpeg")
    # now run the ones you want

    eulerian_view()
    eulerian_vel_quiver_view()
    eulerian_vel_quiver_view(True)
    lagrangian_and_eularian_view()
    velocity_field_view()
    lagrangian_strain_field_view()
    strain_field_view()

demo("toy.jpg")
