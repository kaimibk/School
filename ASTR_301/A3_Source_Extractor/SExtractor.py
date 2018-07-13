from sewpy import SEW
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import ZScaleInterval as ZScale
from shutil import move

class SE():

    def __init__(self, main_dir, out_dir, obj_name, filter):
        print(main_dir)
        self.main_dir = main_dir
        self.out_dir = out_dir
        self.obj_name = obj_name
        self.filter = filter

    def source(self, MINAREA=5, THRESH=1.5, ANALYSIS_THRESH=1.5):
        self.MINAREA  = MINAREA
        self.THRESH = THRESH
        self.ANALYSIS_THRESH = ANALYSIS_THRESH

        sew = SEW(workdir = self.main_dir,

                    params=[
                            "X_IMAGE",
                            "Y_IMAGE",
                            "FLUX_RADIUS(3)",
                            "FLAGS"],

                    config = {
                            ## == Catalog == ##
                            ##"CATALOG_NAME" : "{0}_{1}.cat.txt".format(obj_name,filter),
                            "CATALOG_TYPE" : "ASCII_HEAD",
                            "PARAMETERS_NAME": self.main_dir+"default.param",

                            ## == Extraction == ##
                            "DETECT_TYPE" : "CCD",
                            "DETECT_MINAREA" : self.MINAREA,
                            "DETECT_THRESH" : self.THRESH,
                            "ANALYSIS_THRESH" : self.ANALYSIS_THRESH,

                            "FILTER" : "Y",
                            "FILTER_NAME" : self.main_dir+"default.conv",
                            "DEBLEND_NTHRESH" : 32,
                            "DEBLEND_MINCONT" : 0.005,

                            "CLEAN" : "Y",
                            "CLEAN_PARAM" : 1.0,

                            "MASK_TYPE" : "CORRECT",

                            ## == Photometry == ##
                            "PHOT_APERTURES" : "4,6,8",
                            "PHOT_AUTOPARAMS" : "2.5, 3.5",
                            "SATUR_LEVEL" : 50000.0,

                            "MAG_ZEROPOINT" : 0.0,
                            "MAG_GAMMA" : 4.0,
                            "GAIN" : 0.0,
                            "PIXEL_SCALE" : 0.0,

                            ## == Star/Galaxy Separation == ##
                            "SEEING_FWHM" : 1.2,
                            "STARNNW_NAME" : self.main_dir+"default.nnw",

                            ## == Background == ##
                            "BACK_SIZE" : 64,
                            "BACK_FILTERSIZE" : 3,
                            "BACKPHOTO_TYPE" : "GLOBAL",

                            ## == CHECK IMAGE == ##
                            "CHECKIMAGE_TYPE" : "BACKGROUND,OBJECTS,SEGMENTATION,APERTURES",
                            "CHECKIMAGE_NAME": "{0}{1}_{2}_bg.fits,{0}{1}_{2}_obj.fits,{0}{1}_{2}_seg.fits,{0}{1}_{2}_ap.fits".format(self.out_dir,
                                                                                                                    self.obj_name,
                                                                                                                    self.filter),

                            ## == Memory == ##
                            "MEMORY_OBJSTACK" : 2000,
                            "MEMORY_PIXSTACK" : 200000,
                            "MEMORY_BUFSIZE": 1024,

                            ## == Misc. == ##
                            "VERBOSE_TYPE" : "QUIET"
                            }

                    ##configfilepath = main_dir+"default.sex")
                    )

        out = sew("{0}{1}_{2}.fits".format(self.main_dir, self.obj_name, self.filter))["table"].to_pandas()

        move("{0}{1}_{2}.log.txt".format(self.main_dir, self.obj_name, self.filter),
                 "{0}{1}_{2}.log.txt".format(self.out_dir, self.obj_name, self.filter))

        move("{0}{1}_{2}.cat.txt".format(self.main_dir, self.obj_name, self.filter),
                 "{0}{1}_{2}.cat.txt".format(self.out_dir, self.obj_name, self.filter))


        return out


    def view_obj(self):
        plt.close('all')

        hdu_list = fits.open("{}{}_{}_ap.fits".format(self.out_dir,
                                                    self.obj_name,
                                                    self.filter))
        mpl_fig = plt.figure()

        ax = mpl_fig.add_subplot(111)

        data = hdu_list[0].data

        vmin, vmax = ZScale().get_limits(values=data)

        ax.set_title("MINAREA: {}, DETECT_THRESH: {}, ANALYSIS_THRESH: {}".format(
                                                    self.MINAREA,
                                                    self.THRESH,
                                                    self.ANALYSIS_THRESH))

        ax.imshow(data, cmap="gray",
                         vmin=vmin, vmax=vmax)


        plt.savefig("{}{}_{}_{}-{}-{}.png".format(self.out_dir,
                                                  self.obj_name,
                                                  self.filter,
                                                  self.MINAREA,
                                                  self.THRESH,
                                                  self.ANALYSIS_THRESH),
                    bbox_to_inches="tight", dpi=600)
