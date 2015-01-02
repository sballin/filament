if not os.path.isfile('output/G_matrices/G%d.p'%edge_gridpoints):
    print 'Making G%d.p'%edge_gridpoints
    G = [[0.0 for x in xrange(columns)] for x in xrange(rows)]
    for i,sensor in enumerate(sensors):
        for w in range(edge_gridpoints):
            for h in range(edge_gridpoints):
                r = r_start+w*grid_spacing
                z = z_start+h*grid_spacing
                G[i][h*edge_gridpoints+w] = 1e7*greens_function(r, z-sensor.z, sensor.r, sensor.n_r, sensor.n_z) # 1e7 multiplication for matrix stuff
        progress = i/float(rows)*100
        sys.stdout.write('\r%.2f%%' %progress)
        sys.stdout.flush()
    G_array = numpy.array(G)
    pickle.dump(G_array, open('output/G_matrices/G%d.p'%edge_gridpoints, 'wb'))
    print
else:
    G_array = pickle.load(open('output/G_matrices/G%d.p'%edge_gridpoints, 'rb'))
    print 'DANGER: loaded G%d.p from disk.'%edge_gridpoints
print 'Average value:', numpy.average(G_array)
print 'Matrix rank:', linalg.matrix_rank(G_array)

