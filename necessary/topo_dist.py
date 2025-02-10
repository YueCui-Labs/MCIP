import numpy as np

import nibabel as nib
import scipy.io as scio
from nilearn import surface

# import pygeodesic.geodesic as geodesic
from scipy.sparse import csr_matrix, load_npz, save_npz, diags
import os.path

# use this function to create distance_map.{}.{}.32k_fs_LR.func.gii
def make_topodis(atlas='BN_Atlas'):
    all_vert_num = 32492

    for hemi in ['L', 'R']:

        label_32k = surface.load_surf_data('fsaverage.{}.{}.32k_fs_LR.label.gii'.format(hemi, atlas))
        label_32k[label_32k < 0] = 0
        print(label_32k.shape)

        vertices, faces = surface.load_surf_data('fsaverage.{}.midthickness.32k_fs_LR.surf.gii'.format(hemi))
        print(vertices.shape, faces.shape)
        geoalg = geodesic.PyGeodesicAlgorithmExact(vertices, faces)

        # Create a sparse adjacency matrix
        adj_file = "adjacency_matrix_{}.npz".format(hemi)
        if not os.path.exists(adj_file):
            N = len(vertices)
            adjacency_matrix = np.zeros((N,N), dtype=bool)
            for face in faces:
                # Set the entries of the adjacency matrix based on the distances between connected vertices
                for i in range(3):
                    for j in range(i + 1, 3):
                        vertex1 = face[i]
                        vertex2 = face[j]
                        adjacency_matrix[vertex1, vertex2] = True
                        adjacency_matrix[vertex2, vertex1] = True
            adjacency_matrix = csr_matrix(adjacency_matrix)
            save_npz(adj_file, adjacency_matrix)
        else:
            adjacency_matrix = load_npz(adj_file)

        dists = []
        for ind in np.unique(label_32k[label_32k > 0]):
            ind_vert = np.where(label_32k == ind)[0]
            # print(ind_vert)
            adj = adjacency_matrix[ind_vert,:]
            # print(adj.shape)

            border = []
            for i in range(adj.shape[0]):
                # print(adj[i,:].nonzero()[1])
                vii = label_32k[adj[i,:].nonzero()[1]] != ind  # (1,1,0,0, 0) 1 - is border
                # print(vii)
                if vii.any():
                    border.append(i)

            border = np.array(border, dtype=int)
            # print(border)
            # print(ind_vert[border])

            dist_map, _ = geoalg.geodesicDistances(ind_vert[border], np.arange(all_vert_num, dtype=int))
            print(dist_map.shape)
            dist_map[ind_vert] = - dist_map[ind_vert]
            centdist = np.abs(dist_map.min())
            if centdist == 0:
                centdist = 1
            dist_map /= centdist
            print(dist_map.shape, dist_map.min(), dist_map.max())

            dists.append(dist_map)

        dists = np.stack(dists, axis=1)  # (V, K)
        print(dists.shape)

        template = nib.load("fsaverage.{}.{}.32k_fs_LR.label.gii".format(hemi, atlas))
        template.remove_gifti_data_array(0)
        template.add_gifti_data_array(nib.gifti.gifti.GiftiDataArray(dists.astype(np.float32)))
        nib.loadsave.save(template, 'distance_map.{}.{}.32k_fs_LR.func.gii'.format(hemi, atlas))


def avg_topodis(atlas='BN_Atlas'):
    all_vert_num = 10242

    avg_topodis = []
    for hemi in ['L', 'R']:

        label_32k = surface.load_surf_data('fsaverage.{}.{}.32k_fs_LR.label.gii'.format(hemi, atlas))
        label_32k[label_32k < 0] = 0
        print(label_32k.shape)

        vertices, faces = surface.load_surf_data('fsaverage.{}.midthickness.32k_fs_LR.surf.gii'.format(hemi))
        print(vertices.shape, faces.shape)

        for face in faces:
            # Set the entries of the adjacency matrix based on the distances between connected vertices
            for i in range(3):
                for j in range(i + 1, 3):
                    vertex1 = face[i]
                    vertex2 = face[j]
                    if label_32k[vertex1] >0 and label_32k[vertex2] > 0:
                        avg_topodis.append(np.sqrt(np.sum((vertices[vertex1,:] - vertices[vertex2,:])**2)))

    avg_topodis = np.array(avg_topodis)
    print(avg_topodis.mean(), avg_topodis.std())


if __name__ == '__main__':
    # make_topodis('BN_Atlas')
    make_topodis('Glasser')
    # avg_topodis()
