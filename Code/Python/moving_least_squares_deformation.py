from jax.config import config
config.update("jax_enable_x64", True)
import jax.lax as jl
import jax.numpy as np
from jax import vmap, jacfwd
from bisect import bisect

def cubic_spline(xp, yp):
   """
   cubic spline (not-a-node)
   https://math.libretexts.org/Bookshelves/Applied_Mathematics/Numerical_Methods_(Chasnov)/05%3A_Interpolation/5.03%3A_Cubic_Spline_Interpolation

   Parameters
   ----------
   xp: ndarray with size [nd]
   yp: ndarray with size [nd]

   Assumptions
   -----------
   xp is sorted ascending

   Return
   ------
   annonomous function to evaluate spline
   """

   d = yp[:-1]

   h = xp[1:] - xp[:-1]
   f = yp[1:] - yp[:-1]

   q = f[1:]/h[1:] - f[:-1]/h[:-1]
   q = np.concatenate((np.array([0]), q))
   q = np.concatenate((q, np.array([0])))

   A1 = np.concatenate((h[:-1]/3.0, np.array([-h[-2]-h[-1]])))
   A2 = np.concatenate((np.array([h[1]]), 2.0*(h[:-1]+h[1:])/3.0))
   A2 = np.concatenate((A2, np.array([h[-2]])))
   A3 = np.concatenate((np.array([-h[0]-h[1]]), h[1:]/3.0))
   A4 = np.concatenate((np.array([h[0]]), np.zeros(xp.shape[0]-3)))
   A5 = np.concatenate((np.zeros(xp.shape[0]-3), np.array([h[-1]])))

   A = np.diag(A5, -2) + np.diag(A1,-1) + np.diag(A2,0) + np.diag(A3,1) + np.diag(A4,2)
   b = np.linalg.inv(A) @ q
   c = f/h - (h/3.0)*(b[1:] + 2*b[:-1])
   a = (b[1:] - b[:-1])/(3.0*h)

   b = b[:-1]

   return np.stack((d,c,b,a))

cubic_spline_vec = vmap(cubic_spline, in_axes=(None,0))

def motion_2d_splines(tp, deform):
    """
    compute spline coefficients for 2d motion

    Parameters
    ----------
    tp: ndarray with size [nt]
    deform: ndarray with size [nt, ncp, dim]

    Return
    ------
    deform_coeffs: ndarray with size [ncp*dim, 4, nt]
    """

    nt, ncp, dim = deform.shape
    # rearrange deform to (dim, ncp, nt)
    x_deform = np.moveaxis(deform, (0,1,2),(2,1,0))
    # reshape to (dim*ncp, nt)
    xx_deform = x_deform.reshape(dim*ncp,nt)
    #
    deform_coeffs = cubic_spline_vec(tp, xx_deform)

    return deform_coeffs


def evaluate_spline_one(t, tp, coeffs):
    index = bisect(tp, t) - 1
    if index < 0: index = 0
    if index == len(tp) - 1: index = index - 1
    dt = t - tp[index]
    return coeffs[3,index]*(dt**3) + coeffs[2,index]*(dt**2) + coeffs[1,index]*dt + coeffs[0,index]

evaluate_spline_vec = vmap(evaluate_spline_one, in_axes=(None,None,0))


def mlsd_onepoint(X,ps,qs):
    """
    Moving Least Squares Deformation

    Prameters
    ---------
    X: ndarray with size [dim] reference coordinates for evaluation point
    ps: ndarray with size [n, dim] reference control points
    qs: ndarray with size [n, dim] deformed control points

    Return
    ------
    Deformed x coordinates of evaluation point
    """
    dim = X.shape[0]
    alpha = 1.0
    eps = 1.0e-5
    ws = 1.0 / ((np.linalg.norm(ps - X,axis=1)+eps)**(2.0*alpha))
    wssum = ws.sum()
    pstar = np.dot(ws,ps)/wssum
    qstar = np.dot(ws,qs)/wssum
    phat = ps - pstar
    qhat = qs - qstar

    A = np.zeros((dim,dim))
    pq_outer_products = np.zeros((dim,dim))

    for idx in range(phat.shape[0]):
        A = A + ws[idx]*np.outer(phat[idx,:],phat[idx,:])
        pq_outer_products = pq_outer_products + ws[idx]*np.outer(phat[idx,:],qhat[idx,:])

    M = np.linalg.inv(A) @ pq_outer_products
    return M.T @ (X-pstar) + qstar

mlsd_vec = vmap(mlsd_onepoint, in_axes=(0,None,None), out_axes=0)

mlsdDefGrad = jacfwd(mlsd_onepoint, argnums=0)
mlsdDefGrad_vec = vmap(mlsdDefGrad, in_axes=(0,None,None))

def mlsd_inv_onepoint(x, ps, qs):
   """
   Inverse of Moving Least Squares

   Computed using reversed MLS for guess and then nonlinear solver
   """

   def f(X):
      return x - mlsd_onepoint(X, ps, qs)
   def df(X):
      return -mlsdDefGrad(X, ps, qs)
   def fsolve(f, X0, df):
      def cond(X):
         return np.linalg.norm(f(X)) > 1.49012e-08
      def body(X):
         return X - np.linalg.inv(df(X)) @ f(X)

      X = jl.while_loop(cond, body, X0)
      return X

   # use reversed control points with mlsd to generate guess for solution
   X_guess = mlsd_onepoint(x, qs, ps)
   # polish guess to get good solution
   X = fsolve(f, X_guess, df)
   return X

mlsd_inv_vec = vmap(mlsd_inv_onepoint, in_axes=(0,None,None), out_axes=0)


def current_2d_deform(t, tp, deform_coeffs):
   dim = 2
   vv = evaluate_spline_vec(t, tp, deform_coeffs)
   vv = np.swapaxes(vv.reshape(dim,-1), 0,1)
   return vv

def mlsd_2d_motion_onepoint(X, t, tp, deform_coeffs):
   """
   Motion with MLS

   Parameters
   ----------
   X: ndarray with size [dim] reference coordinates for evaluation point
   t: evaluation time
   tp: ndarray with size [nt] control times for cubic splines
   deform_coeffs: ndarray with size [n*dim, 4, nt] of cubic spline coeffs for deformation

   Return
   ------
   Deformed x coordinates at evaluation point,time
   """

   deform_ref = current_2d_deform(0.0, tp, deform_coeffs)
   deform = current_2d_deform(t, tp, deform_coeffs)

   x = mlsd_onepoint(X, deform_ref, deform)
   return x

def mlsd_2d_inv_motion_onepoint(x, t, tp, deform_coeffs):
   """
   Inversion Motion with MLS

   Parameters
   ----------
   x: ndarray with size [dim] current coordinates for evaluation point
   t: evaluation time
   tp: ndarray with size [nt] control times for cubic splines
   deform_coeffs: ndarray with size [n*dim, 4, nt] of cubic spline coeffs for deformation

   Return
   ------
   Reference X coordinates at evaluation point,time
   """

   deform_ref = current_2d_deform(0.0, tp, deform_coeffs)
   deform = current_2d_deform(t, tp, deform_coeffs)

   X = mlsd_inv_onepoint(x, deform_ref, deform)
   return X

mlsd_2d_motion_defgrad_onepoint = jacfwd(mlsd_2d_motion_onepoint, argnums=0)
mlsd_2d_motion_defgrad_dot_onepoint = jacfwd(mlsd_2d_motion_defgrad_onepoint, argnums=1)
mlsd_2d_motion_dot_onepoint = jacfwd(mlsd_2d_motion_onepoint, argnums=1)
def mlsd_2d_motion_vel_onepoint(x, t, tp, deform_coeffs):
   X = mlsd_2d_inv_motion_onepoint(x, t, tp, deform_coeffs)
   return mlsd_2d_motion_dot_onepoint(X, t, tp, deform_coeffs)
mlsd_2d_motion_velgrad_onepoint = jacfwd(mlsd_2d_motion_vel_onepoint, argnums=0)

mlsd_2d_motion_vec = vmap(mlsd_2d_motion_onepoint, in_axes=(0,None,None,None), out_axes=0)
mlsd_2d_inv_motion_vec = vmap(mlsd_2d_motion_onepoint, in_axes=(0,None,None,None), out_axes=0)
mlsd_2d_motion_defgrad_vec = vmap(mlsd_2d_motion_defgrad_onepoint, in_axes=(0,None,None,None), out_axes=0)
mlsd_2d_motion_defgrad_dot_vec = vmap(mlsd_2d_motion_defgrad_dot_onepoint, in_axes=(0,None,None,None), out_axes=0)
mlsd_2d_motion_dot_vec = vmap(mlsd_2d_motion_dot_onepoint, in_axes=(0,None,None,None), out_axes=0)
mlsd_2d_motion_vel_vec = vmap(mlsd_2d_motion_vel_onepoint, in_axes=(0,None,None,None), out_axes=0)
mlsd_2d_motion_velgrad_vec = vmap(mlsd_2d_motion_velgrad_onepoint, in_axes=(0,None,None,None), out_axes=0)
