def nilearn_subplot(fig_params, plot_params):
    """
    Create multiview whole-brain surface subplots with a common colorbar displayed
    """
    
    # Create global colorbar
    cmap = plot_params['cmap']
    data = plot_params['data']
    threshold = plot_params['threshold']
    
    # Create subplots
    dims = fig_params['dim']
    figsize = fig_params['figsize']
    title = fig_params['title']
    fname = fig_params['fname']
    surf_mesh = plot_params['surf_mesh']
    bg_map = plot_params['bg_map']
    hemis = plot_params['hemis']
    views = plot_params['views']
    symmetric = plot_params['symmetric']
    
    fig, axes = plt.subplots(nrows=dims[0], ncols=dims[1], subplot_kw={'projection': '3d'}, figsize=figsize)
                             
    ax_flat = axes.flat

    vmin = 1
    vmax = 4
    
    our_cmap = get_cmap(cmap)
    norm = Normalize(vmin=vmin, vmax=vmax)
#     if threshold is not None:  # threshold current cmap
#         cmaplist = [our_cmap(i) for i in range(our_cmap.N)]
#         istart = int(norm(-threshold, clip=True) * (our_cmap.N - 1))
#         istop = int(norm(threshold, clip=True) * (our_cmap.N - 1))
#         for i in range(istart, istop):
#             cmaplist[i] = (0.5, 0.5, 0.5, 1.)
#         our_cmap = LinearSegmentedColormap.from_list('Custom cmap', cmaplist, our_cmap.N)

    for i, ax in enumerate(ax_flat):
        nilearn.plotting.plot_surf(
            surf_mesh=surf_mesh[i],
            surf_map=data[:int(len(data) / 2)] if hemis[i] == 'left' else data[int(len(data) / 2):],
            hemi=hemis[i],
            view=views[i], 
            bg_map=bg_map[i], cbar_vmin=vmin, cbar_vmax=vmax, cmap=cmap, 
            symmetric_cmap=symmetric, vmin=vmin, vmax=vmax, threshold=threshold, axes=ax
        )
                             
        ax.dist = 6  # Decrease to zoom-in

    proxy_mappable = ScalarMappable(cmap=our_cmap, norm=norm)
    proxy_mappable.set_array(data)

    nb_ticks = 5
    ticks = np.linspace(vmin, vmax, nb_ticks)
    bounds = np.linspace(vmin, vmax, our_cmap.N)
    
    p0 = ax_flat[len(ax_flat) - dims[1]].get_position().get_points().flatten()
    p1 = ax_flat[len(ax_flat) - 1].get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], 0, p1[2]-p0[0], 0.025])
    plt.colorbar(proxy_mappable, cax=ax_cbar, boundaries=bounds, ticks=ticks, spacing='proportional', orientation='horizontal', format='%.0e')
    
    fig.text(0.5, 1, title, ha='center', va='center', fontsize=34)
    
    plt.tight_layout()
    if fname:
        plt.savefig(fname, bbox_inches='tight')