import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        return np.min(zs)

def generate_cylinder(r, h, m_twist, tilt_angle, res=50):
    z = np.linspace(0, h, res)
    theta = np.linspace(0, 2*np.pi, res)
    theta_grid, z_grid = np.meshgrid(theta, z)
    
    # Twist factor
    twist = m_twist * (z_grid / h) * 2 * np.pi
    
    x_grid = r * np.cos(theta_grid + twist)
    y_grid = r * np.sin(theta_grid + twist)
    
    # Rotation (Tilt around Y-axis)
    rad = np.radians(tilt_angle)
    c, s = np.cos(rad), np.sin(rad)
    
    x_rot = x_grid * c + z_grid * s
    y_rot = y_grid
    z_rot = -x_grid * s + z_grid * c
    
    return x_rot, y_rot, z_rot

def plot_fig1():
    fig = plt.figure(figsize=(12, 6))
    
    # Parameters
    R = 1.0
    H = 4.0
    M_twist = 1.0 # Visual twist
    Tilt = 15.0   # Degrees
    
    # --- Left Panel: Aligned ---
    ax1 = fig.add_subplot(121, projection='3d')
    x, y, z = generate_cylinder(R, H, M_twist, 0)
    
    # Bulk Planes (Top and Bottom)
    xx, yy = np.meshgrid(np.linspace(-2, 2, 10), np.linspace(-2, 2, 10))
    zz_bottom = np.zeros_like(xx)
    zz_top = np.full_like(xx, H)
    
    # Plot Bulk
    ax1.plot_surface(xx, yy, zz_bottom, alpha=0.2, color='gray')
    ax1.plot_surface(xx, yy, zz_top, alpha=0.2, color='gray')
    
    # Plot Handle
    ax1.plot_surface(x, y, z, alpha=0.6, color='cyan', edgecolor='none')
    
    # Twist lines
    for i in range(0, 50, 5):
        ax1.plot(x[:,i], y[:,i], z[:,i], color='blue', linewidth=0.5, alpha=0.5)

    # Axis Arrow
    a = Arrow3D([0, 0], [0, 0], [-0.5, H+0.5], mutation_scale=20, arrowstyle="-|>", color="k", lw=2)
    ax1.add_artist(a)
    
    # Labels
    ax1.text(0, 0, H+1, r"Axis $e^1$", ha='center')
    ax1.set_title(r"(a) Aligned ($\varepsilon=0$)", fontsize=14, y=1.05)
    ax1.set_axis_off()
    ax1.set_box_aspect([1,1,2])
    
    # --- Right Panel: Precessed ---
    ax2 = fig.add_subplot(122, projection='3d')
    x_t, y_t, z_t = generate_cylinder(R, H, M_twist, Tilt)
    
    # Plot Bulk (Fixed geometry)
    ax2.plot_surface(xx, yy, zz_bottom, alpha=0.2, color='gray')
    ax2.plot_surface(xx, yy, zz_top, alpha=0.2, color='gray')
    
    # Plot Handle (Tilted)
    ax2.plot_surface(x_t, y_t, z_t, alpha=0.6, color='magenta', edgecolor='none')
    
    # Twist lines
    for i in range(0, 50, 5):
        ax2.plot(x_t[:,i], y_t[:,i], z_t[:,i], color='darkred', linewidth=0.5, alpha=0.5)

    # Axis Arrow (Tilted)
    # Calculate tip position
    rad = np.radians(Tilt)
    tip_x = (H+0.5) * np.sin(rad)
    tip_z = (H+0.5) * np.cos(rad)
    
    a2 = Arrow3D([0, tip_x], [0, 0], [0, tip_z], mutation_scale=20, arrowstyle="-|>", color="k", lw=2)
    ax2.add_artist(a2)
    
    # Mismatch Annotation
    # Circle on Bulk (Expected) vs Ellipse (Actual)
    theta_circ = np.linspace(0, 2*np.pi, 100)
    xc = R * np.cos(theta_circ)
    yc = R * np.sin(theta_circ)
    ax2.plot(xc, yc, np.full_like(xc, H), 'k--', linewidth=1.5, label='Bulk Geometry')
    
    # Top rim of tilted cylinder
    x_rim = x_t[-1,:]
    y_rim = y_t[-1,:]
    z_rim = z_t[-1,:]
    ax2.plot(x_rim, y_rim, z_rim, 'r-', linewidth=2, label='Handle Geometry')

    # Epsilon Angle
    ax2.text(0, 0, H+1, r"Tilt $\varepsilon(t)$", ha='center')
    ax2.text(1.5, 1.5, H, r"Boundary Mismatch $\to \Lambda_q \varepsilon$", color='red', fontsize=10)

    ax2.set_title(r"(b) Precessed ($\varepsilon \neq 0$)", fontsize=14, y=1.05)
    ax2.set_axis_off()
    ax2.set_box_aspect([1,1,2])

    # Save
    plt.tight_layout()
    plt.savefig('Fig1_Junction_Mismatch.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    plot_fig1()

