#!/usr/bin/env python3
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from scipy.spatial.distance import cdist
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from optparse import OptionParser
#plt.tight_layout()

def smear(data, sigma):
    """
    Apply Gaussian smearing to spectrum y value.

    Args:
        sigma: Std dev for Gaussian smear function
    """
    diff = [data[0, i + 1] - data[0, i] for i in range(np.shape(data)[0] - 1)]
    avg_x_per_step = np.sum(diff) / len(diff)
    data[1, :] = gaussian_filter1d(data[1, :], sigma / avg_x_per_step)
    return data


class RDF(object):
    """a class of crystal structure.
    Args:
        crystal: crystal class from pymatgen
        symmetrize: symmetrize the structure before RDF computation
        R_max: maximum cutoff distance
        R_bin: length of each bin when computing the RDF
        width: width of gaussian smearing
    Attributes:
        crystal
        R_max
        R_bin
        width
        RDF
        plot_RDF()
    """

    def __init__(self, crystal, symmetrize=True,
                 R_max=6, R_bin=0.2, sigma=0.2):
        """Return a RDF object with the proper info"""
        self.R_max = R_max
        self.R_bin = R_bin
        self.sigma = sigma
        if symmetrize:
            finder = SpacegroupAnalyzer(crystal, symprec=0.06,
                                        angle_tolerance=5)
            crystal = finder.get_conventional_standard_structure()

        self.compute_RDF(crystal)
        # self.plot_RDF()

    def compute_RDF(self, crystal):
        """
        Computes the radial distribution function of a given crystal.
        Args:
        self: RDF
        crystal: Crystal structure information
        Returns: None
        """
        R_max = self.R_max
        R_bin = self.R_bin

        # below is the old code before vectorization

        # R: the array which contains the number of occurences of atomic pairs
        # in each [dmin, dmax].

        # rij_dot: the atomic coordinates of supercell in cartesian format

        # rij_dist: the distance matrix between atoms in the big cell and in
        # the small cell
        # the idea is
        # 1, to loop over all atoms in the small cell
        # 2, calculate rij_dist
        # 3, for each distance bin [dmin, dmax], count the occurences of
        # distances

        # R = np.zeros(round(R_max/R_bin))
        # rij_dot = self.find_supercell(crystal, R_max)
        # for atom in crystal.frac_coords:
        #    origin = [np.dot(atom, crystal.lattice.matrix)]
        #    rij_dist = cdist(rij_dot, origin)
        #    for i in range(len(R)):
        #        d_min, d_max = R_bin*(i+0.5), R_bin*(i+1.5)
        #        R[i] += len([x for x in rij_dist if d_min<=x<d_max])

        # vectotrized version:
        # Vectorizing code refers to operations that are performed
        # on multiple components of a vector from a single statement

        # apply_along_axis applies a function across the dimension of an array
        # in essence, it is an optimized for loop.

        # see find supercell method
        rij_dot = self.find_supercell(crystal, R_max) # 原来晶胞的的
        length = round(R_max/R_bin)  # length of distance array
        # create minimum distance vector and add dimension
        d_min = np.arange(0.5, length+0.5, 1)
        # create maximum distance vector and add dimension
        d_max = np.arange(1.5, length+1.5, 1)
        # stack the min and max distance vectors into array
        d = np.vstack((d_min, d_max))*R_bin

        def compute_rij_dist(atom):
            """
            Computes the distance between atoms in the unit cell and atoms in
            the supercell
            Args:
            atom = Fractional coordinates of the atom within the unit cell
            Returns: the euclidean distance between atoms in the unit cell and
                     atoms in the supercell
            """
            # dot product of atomic fractional coordinates and lattice matrix
            # print(atom); input("atom")
            origin = np.dot(atom, crystal.lattice.matrix) # atom = [0.  0.  0.5]  origin = [0.       0.       2.820871]
            # print(origin); input("origin")
            origin = origin[np.newaxis, :]  # add dimension to array: origin = [[0.       0.       2.820871]]
            # print(origin); input("origin np.newaxis")
            # print(rij_dot, rij_dot.shape); input("rij_dot")
            # x = cdist(rij_dot, origin)
            # print(x, len(x)); input("cdist(rij_dot, origin)")
            return cdist(rij_dot, origin)
            # return x
        
        # loop over fractional coordinates of each atom in the crystal to
        # compute an array of euclidean distances
        # print(crystal.frac_coords); input("crystal.frac_coords")
        rij_dist = np.apply_along_axis(compute_rij_dist, axis=1, # 沿着第2个轴（列）应用函数
                                       arr=crystal.frac_coords)
        # print(rij_dist.shape); input("rij_dist")

        def compute_R(span):
            """
            Counts the euclidean distances within a bin range
            Args:
            span = An ordered pair of min and max distances for the bin
            Returns: an count of distances within the bin range
            """
            return ((span[0] <= rij_dist) & (rij_dist < span[1])).sum()

        # R: the array which contains the number of occurences of atomic pairs
        # in each [dmin, dmax].
        R = np.apply_along_axis(compute_R, axis=0, arr=d)

        # radii in angstrom
        r = np.arange(1, length+1, 1)*R_bin

        # now calculate RDF based on the equation *** (reference from the book)
        r = np.arange(1, length+1, 1)*R_bin
        rho = len(crystal.frac_coords)/crystal.volume
        R = R/(4*np.pi*R_bin*rho**2*crystal.volume) * np.power(r, -2)
        self.RDF = np.vstack((r, R))
        return

    def plot_RDF(self, filename=None, tag=None):
        """ plot PXRD """
        plt.rcParams["font.weight"] = "bold"
        plt.rcParams["axes.labelweight"] = "bold"
        plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
        plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
        datax = smear(self.RDF, self.sigma)
        #print("self.RDF",self.RDF)
        #plt.plot(datax[0, :], datax[1, :])
        plt.plot(self.RDF[0, :],self.RDF[1, :],label=tag,linewidth=2.0)
        plt.xlabel("$r(\AA)$",fontsize=12)
        plt.ylabel("$g(r)$",fontsize=12)
        if filename is None:
            plt.show()
            #pass
        else:
            plt.savefig(filename)
            plt.close()

    def find_supercell(self, crystal, R_max, atom=[0.5, 0.5, 0.5]):
        """
        Finds the supercell in cartesian format
        Args:
        crystal: Structure
        R_max: maximum radii of supercell
        atom: index size of supercell
        Returns: the atomic coordinates of supercell in cartesian format
        """

        def calculateR(vect):
            return np.linalg.norm(np.dot(vect-atom, crystal.lattice.matrix))

        # temporary max index
        hkl_max = np.array([1, 1, 1])

        # cartesian product of [-1,0,1] with itself returns all possible index
        # configurations with no repeats
        hkl_index = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3)

        # calculate R value across all rows of an array
        R = np.apply_along_axis(calculateR, axis=1, arr=hkl_index)

        # Creates multiple array from R
        multiple = np.round(R_max/R).astype(dtype='float64')

        # multiplies the each row in the index by the corresponding multiple
        hkl_index = (hkl_index.T * multiple).T

        # finds the max indeces
        hkl_max[0] = np.amax(hkl_index[:, 0])
        hkl_max[1] = np.amax(hkl_index[:, 1])
        hkl_max[2] = np.amax(hkl_index[:, 2])
        h, k, l = hkl_max

        # creates arrays for the integers between the max indeces
        i = np.arange(-h, h+1)
        j = np.arange(-k, k+1)
        k = np.arange(-l, l+1)

        # cartesian product or arrays i, j and k
        ijk_index = np.array(np.meshgrid(i, j, k)).T.reshape(-1, 3)

        # shapes of crystal fractional coordinates and ijk index array
        fracCoordsSize = np.shape(crystal.frac_coords)
        ijk_Size = np.shape(ijk_index)

        coors = np.empty([fracCoordsSize[0]*ijk_Size[0], 3])

        for c, coor in enumerate(crystal.frac_coords):
            coors[c*ijk_Size[0]:c*ijk_Size[0]+ijk_Size[0]] = ijk_index+coor

        supercell = np.asarray(coors)
        rij_dot = np.dot(supercell, crystal.lattice.matrix)
        return rij_dot


if __name__ == "__main__":
    # -------------------------------- Options -------------------------
    parser = OptionParser()
    parser.add_option("-c", "--crystal", dest="structure", default='',
                      help="crystal from file, cif or poscar, REQUIRED")
    parser.add_option("-r", "--Rmax", dest="Rmax", default=6, type='float',
                      help="Rmax, default: 6 A", metavar="Rmax")
    parser.add_option("-s", "--sigma", dest="sigma", default=0.20, type='float',
                      help="sigma width for Gaussian smearing, default: 0.20")
    parser.add_option("-d", "--delta", dest="delta", default=0.08,
                      type='float', help="step length, default: 0.08")
    parser.add_option("-o", "--output", dest="mstyle", default='bmh',
                      help="matplotlib style, fivethirtyeight, bmh, grayscale, dark_background, ggplot")
    parser.add_option("-p", "--plot", dest="plot", default=None,
                      help="generate the plot to file, default: None")

    (options, args) = parser.parse_args()
    #plt.style.use('classic')
    fig=plt.figure(figsize=(8,6))
    
    #ax = plt.axes()
    #ax.set_facecolor("white")
    test = Structure.from_file(options.structure)
    infile=os.popen('cat %s' %(options.structure))
    tag=infile.readline()[:-1]
    rdf = RDF(test, R_max=options.Rmax, R_bin=options.delta,
                  sigma=options.sigma)
    #print('-----RDF value-----')
    #print(rdf.RDF)
    rdf.plot_RDF(options.plot,tag=tag)
    #lt.rcParams.update({'font.size': 20})#设置图例大小
    plt.xticks(size=12)#刻度线大小
    plt.yticks(size=12)
    plt.legend(bbox_to_anchor=(0.80,0.40),#图例边界框起始位置
               fontsize=12,#设置图例大小
               #loc="lower right",#图例位置
               ncol=1,#列数
               mode="None",#当前位置为"expend"时，图例会扩展到整个坐标轴区域
               borderaxespad=0,#坐标轴和图例边界的间距
               #title="method",#图例标题
               shadow=False,#是否为线框添加阴影
               fancybox=True)#线框圆角处理参数
    plt.show()
    plt.savefig('./gor.png')

