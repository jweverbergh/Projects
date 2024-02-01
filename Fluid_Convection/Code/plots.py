import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter


def anim(prefix, n_start=0, save=False):
    
    def update(n_idx):
        with open("./save_txt/" +prefix + str(n_idx+n_start) + ".txt", "r") as file:
            Theta = np.loadtxt(file, max_rows=ny) / DT
            #w = np.loadtxt(file) * H / U
            u = np.loadtxt(file, max_rows=ny) / U
            v = np.loadtxt(file, max_rows=ny+1) / U
            w = ((v[1:ny, 1:] - v[1:ny, :-1]) - (u[1:, 1:nx] - u[:-1, 1:nx])) / h * H
            norm_u = 0.5 * np.sqrt((u[:, 1:] + u[:, :-1])**2 + (v[1:, :] + v[:-1, :])**2)
        time_text.set_text(time_template.format((n_idx+n_start) / frame_save))
        pcolormsh_T.set_array(Theta)
        pcolormsh_w.set_array(w)
        pcolormsh_norm_u.set_array(norm_u)
        return time_text, pcolormsh_T, pcolormsh_w


    with open("./save_txt/" + prefix + "_param.txt", "r") as file:
        L, H, DT, U, nu, beta, qw, k, g, alpha, l0, Gr, Pr, nx, ny, h, dt, nt, frame_save, n_save = [el for el in np.loadtxt(file, max_rows=1)]

    nx = int(nx); ny = int(ny); nt = int(nt); n_save = int(n_save);
    print(nx, ny)
    x_theta = np.linspace(0, L, nx)
    y_theta = np.linspace(0, H, ny)
    x_w = np.linspace(0, L-h, nx-1)
    y_w = np.linspace(0, H-h, ny-1)
    print("dt = {}, t_max = {}".format(dt * U / H, nt * dt * U / H))
    
    fig, axs = plt.subplots(1, 3, figsize=(14., 9.))
    with open("./save_txt/" +prefix + "maxVals.txt", "r") as file:
        Theta_MAX, w_MAX, w_min, normu_MAX = [el for el in np.loadtxt(file)]

    pcolormsh_T = axs[0].pcolormesh(x_theta, y_theta, np.zeros((ny, nx)), cmap=plt.get_cmap('hot'), vmin=0, vmax=Theta_MAX)#, shading='gouraud') #, **kw)
    plt.colorbar(pcolormsh_T, ax=axs[0], orientation='horizontal') 
    pcolormsh_w = axs[1].pcolormesh(x_w, y_w, np.zeros((ny-1, nx-1)), cmap=plt.set_cmap('bwr'), vmin=-1, vmax=1)#, shading='gouraud') #, **kw)
    plt.colorbar(pcolormsh_w, ax=axs[1], orientation='horizontal') 
    pcolormsh_norm_u = axs[2].pcolormesh(x_theta, y_theta, np.zeros((ny, nx)), cmap=plt.set_cmap('jet'), vmin=0, vmax=normu_MAX)#, shading='gouraud') #, **kw)
    plt.colorbar(pcolormsh_norm_u, ax=axs[2], orientation='horizontal')

    time_template = r'$t = \mathtt{{{:.2f}}} \;[s]$'
    time_text = plt.gcf().text(0.01, 0.975, 'test', fontsize=13)
    axs[0].set_aspect("equal")
    axs[1].set_aspect("equal")
    axs[2].set_aspect("equal")
    axs[0].set_title("$\\frac{T-T_{\\infty}}{\\Delta T}$")
    axs[1].set_title("$\\frac{\\omega H}{U}$")
    axs[2].set_title("$\\frac{|v|}{U}$")
    fig.tight_layout()

    anim = FuncAnimation(fig, update, n_save-n_start, init_func=lambda: None, interval=1/frame_save, repeat_delay=2000)

    if save:
        plt.show()
        """         writerMP4 = FFMpegWriter(fps=24)
        anim.save("./Animation/anim_" + prefix + ".mp4", writer=writerMP4) """
        writergif = PillowWriter(fps=24)
        anim.save("./Animation/anim_" + prefix + ".gif", writer=writergif)
    else:
        plt.show()



def plot_fields(prefix, idx):
    with open("./save_txt/" + prefix + "_param.txt", "r") as file:
        L, H, DT, U, nu, beta, qw, k, g, alpha, l0, Gr, Pr, nx, ny, h, dt, nt, frame_save, n_save = [el for el in np.loadtxt(file, max_rows=1)]

    nx = int(nx); ny = int(ny); nt = int(nt); n_save = int(n_save)
    print(nx, ny)
    x_theta = np.linspace(0, L, nx)
    y_theta = np.linspace(0, H, ny)
    x_w = np.linspace(0, L-h, nx-1)
    y_w = np.linspace(0, H-h, ny-1)
    print("dt = {}, t_max = {}".format(dt * U / H, nt * dt * U / H))

    fig, axs = plt.subplots(1, 3, figsize=(14., 9.))
    with open("./save_txt/" +prefix + "maxVals.txt", "r") as file:
        Theta_MAX, w_MAX, w_min, normu_MAX = [el for el in np.loadtxt(file)]


    pcolormsh_T = axs[0].pcolormesh(x_theta, y_theta, np.zeros((ny, nx)), cmap=plt.get_cmap('hot'), vmin=0, vmax=Theta_MAX, shading='gouraud') #, **kw)
    cbar1 = plt.colorbar(pcolormsh_T, ax=axs[0], orientation='horizontal') 
    cbar1.ax.tick_params(labelsize=16)
    pcolormsh_w = axs[1].pcolormesh(x_w, y_w, np.zeros((ny-1, nx-1)), cmap=plt.set_cmap('bwr'), vmin=-1, vmax=1, shading='gouraud') #, **kw)
    cbar2 = plt.colorbar(pcolormsh_w, ax=axs[1], orientation='horizontal')
    cbar2.ax.tick_params(labelsize=16)
    pcolormsh_norm_u = axs[2].pcolormesh(x_theta, y_theta, np.zeros((ny, nx)), cmap=plt.set_cmap('jet'), vmin=0, vmax=normu_MAX, shading='gouraud') #, **kw)
    cbar3 = plt.colorbar(pcolormsh_norm_u, ax=axs[2], orientation='horizontal')
    cbar3.ax.tick_params(labelsize=16)

    with open("./save_txt/" +prefix + str(idx) + ".txt", "r") as file:
        Theta = np.loadtxt(file, max_rows=ny) / DT
        print(np.shape(Theta))
        u = np.loadtxt(file, max_rows=ny) / U
        v = np.loadtxt(file, max_rows=ny+1) / U
        w = ((v[1:ny, 1:] - v[1:ny, :-1]) - (u[1:, 1:nx] - u[:-1, 1:nx])) / h * H
        print(np.shape(u[:, 1:] + u[:, :-1]))
        print(np.shape(v[1:, :] + v[:-1, :]))
        norm_u = 0.5 * np.sqrt((u[:, 1:] + u[:, :-1])**2 + (v[1:, :] + v[:-1, :])**2)
        print(np.shape(norm_u))
        pcolormsh_T.set_array(Theta)
        pcolormsh_w.set_array(w)
        pcolormsh_norm_u.set_array(norm_u)

    axs[0].set_aspect("equal")
    axs[1].set_aspect("equal")
    axs[2].set_aspect("equal")
    axs[0].set_title("$\\frac{T-T_{\\infty}}{\\Delta T}$", fontsize=22)
    axs[1].set_title("$\\frac{\\omega H}{U}$", fontsize=22)
    axs[2].set_title("$\\frac{|v|}{U}$", fontsize=22)
    axs[0].tick_params(axis='both', labelsize=16)
    axs[1].tick_params(axis='both', labelsize=16)
    axs[2].tick_params(axis='both', labelsize=16)
    fig.tight_layout()
    plt.savefig("./save_plot/" + prefix + "_fields_at_adim_time_" + str(int(idx/2*10)) + ".pdf", format="pdf")
    plt.show()
    


def convert_plot_values(prefix):
    with open("./save_txt/" + prefix + "_param.txt", "r") as file:
        L, H, DT, U, nu, beta, qw, k, g, alpha, l0, Gr, Pr, nx, ny, h, dt, nt, frame_save, n_save = [el for el in np.loadtxt(file, max_rows=1)]

    D = 0.6 * L
    nx = int(nx); ny = int(ny); nt = int(nt); n_save = int(n_save);
    print(nx, ny)
    print("dt = {}, t_max = {}".format(dt * U / H, nt * dt * U / H))
    times = np.arange(n_save+1) / frame_save
    x = np.arange(nx)*h + h/2 - L/2
    y = np.arange(ny)*h + h/2 - L/2

    mask = (np.sqrt(np.outer(x, np.ones(ny)).T**2 + np.outer(np.ones(nx), y).T**2) ) <= D / 2
    # test mask
    fig, ax = plt.subplots(1, 1, figsize=(5, 7.5))
    ax.pcolormesh(x+L/2, y+L/2, np.zeros((ny, nx)), cmap=plt.get_cmap('hot'), vmin=0, vmax=1).set_array(mask)
    plt.show()

    Theta_avg = np.zeros(n_save+1)
    flux_ratio = np.zeros(n_save+1)
    Theta_cyl = np.zeros(n_save+1)
    Theta_rms = np.zeros(n_save+1)
    Re_h = np.zeros(n_save+1)
    Re_w = np.zeros(n_save+1)

    Theta_MAX = 0
    w_MAX = 0
    w_min = 0
    normu_MAX = 0

    for i in range(n_save + 1):
        with open("./save_txt/" +prefix + str(i) + ".txt", "r") as file:
            Theta = np.loadtxt(file, max_rows=ny) / DT
            u = np.loadtxt(file, max_rows=ny) / U
            v = np.loadtxt(file, max_rows=ny+1) / U
            w = ((v[1:ny, 1:] - v[1:ny, :-1]) - (u[1:, 1:nx] - u[:-1, 1:nx])) / h * H
            uv_max = 0.5*np.max(np.abs(u[:, 1:] + u[:, :-1])) + 0.5*np.max(np.abs(v[1:, :] + v[:-1, :]))
            norm_u = 0.5 * np.sqrt((u[:, 1:] + u[:, :-1])**2 + (v[1:, :] + v[:-1, :])**2)
        Theta_avg[i] = np.sum(Theta) * h**2 / (L * H)
        Theta_cyl[i] = np.sum(Theta * mask) * h**2 / (np.pi * (D/2)**2)
        flux_ratio[i] = np.sum(Theta[-1]) * k / (l0 + 0.5*h) * h / L / qw
        Theta_rms[i] = np.sqrt(np.sum((Theta - Theta_avg[i])**2) * h**2 / (L * H))
        Re_h[i] = uv_max * h / nu
        Re_w[i] = np.max(np.abs(w)) * h**2 / nu

        Theta_MAX = max(Theta_MAX, np.max(Theta))
        w_MAX = max(w_MAX, np.max(w))
        w_min = min(w_min, np.min(w))
        normu_MAX = max(normu_MAX, np.max(norm_u))

    
    np.savetxt("./save_txt/" +prefix + "plots.txt", (Re_h, Re_w, flux_ratio, Theta_avg, Theta_cyl, Theta_rms))
    print((Theta_MAX, w_MAX, w_min, normu_MAX))
    np.savetxt("./save_txt/" +prefix + "maxVals.txt", (Theta_MAX, w_MAX, w_min, normu_MAX))


def plot_temperatures(prefix):
    with open("./save_txt/" + prefix + "_param.txt", "r") as file:
        L, H, DT, U, nu, beta, qw, k, g, alpha, l0, Gr, Pr, nx, ny, h, dt, nt, frame_save, n_save = [el for el in np.loadtxt(file, max_rows=1)]
    times = np.arange(n_save+1) / frame_save
    values = np.loadtxt("./save_txt/" + prefix + "plots.txt")

    fig, axs = plt.subplots(2, 1, figsize=(10., 10.))
    axs[0].plot(times, values[0], label='')
    axs[1].plot(times, values[1], label=(''))
    axs[1].set_xlabel("$\\frac{tU}{H}$", fontsize=20)
    axs[0].set_ylabel("$Re_h$", fontsize=18)
    axs[1].set_ylabel('$Re_\\omega$', fontsize=18)
    axs[0].tick_params(axis='both', labelsize=16)
    axs[1].tick_params(axis='both', labelsize=16)
    for i in range(2):
        axs[i].grid(ls=':')
    axs[0].set_title("Mesh Reynolds numbers", fontsize=18)
    fig.tight_layout()
    plt.savefig("./save_plot/" + prefix + "_reynolds.pdf", format="pdf")
    plt.show()
    

    fig, ax = plt.subplots(1, 1, figsize=(10., 5.))
    ax.plot(times, values[2])
    ax.set_xlabel("$\\frac{tU}{H}$", fontsize=20)
    ax.set_ylabel("$\\frac{<q_e>(t)}{q_w}$", fontsize=18)
    ax.set_title("Heat flux density ratio", fontsize=18)
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(ls=':')
    fig.tight_layout()
    plt.savefig("./save_plot/" + prefix + "_heat_flux_density_ratio.pdf", format="pdf")
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(10., 5.))
    ax.plot(times, values[3])
    ax.set_xlabel("$\\frac{tU}{H}$", fontsize=20)
    ax.set_ylabel("$\\frac{<T>(t) - T_{\\infty}}{\\Delta T}$", fontsize=18)
    ax.set_title("Dimensionless spatially-averaged fluid temperature", fontsize=18)
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(ls=':')
    fig.tight_layout()
    plt.savefig("./save_plot/" + prefix + "_fluid_temp.pdf", format="pdf")
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(10., 5.))
    ax.plot(times, values[4])
    ax.set_xlabel("$\\frac{tU}{H}$", fontsize=20)
    ax.set_ylabel("$\\frac{<T>|_{cyl}(t) - T_{\\infty}}{\\Delta T}$", fontsize=18)
    ax.set_title("Dimensionless temperature of the mixer and nearby fluid", fontsize=18)
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(ls=":")
    fig.tight_layout()
    plt.savefig("./save_plot/" + prefix + "_mixer_temp.pdf", format="pdf")
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(10., 5.))
    ax.plot(times, values[5])
    ax.set_xlabel("$\\frac{tU}{H}$", fontsize=20)
    ax.set_ylabel("$\\frac{T_{rms}(t) - T_{\\infty}}{\\Delta T}$", fontsize=18)
    ax.set_title("Dimensionless rms of the fluid temperature fluctuations", fontsize=18)
    ax.tick_params(axis='both', labelsize=16)
    ax.grid(ls=':')
    fig.tight_layout()
    plt.savefig("./save_plot/" + prefix + "_rms_temp.pdf", format="pdf")
    plt.show()


#anim("mixer", n_start=0, save=False)
plot_fields("mixer", 10)
# convert_plot_values("mixer")
# plot_temperatures("mixer")