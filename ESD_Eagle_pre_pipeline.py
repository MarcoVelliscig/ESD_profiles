#/usr/bin/python

# the ESD object contain the ESD profiles (either read from data files or
# computed autonomously)
# the satellite fractions
# all the plotting options
# statistics for various quantities in the stellar mass bins (for simulations)

class ESD(object):
    def __init__(self,  name , mass_arr , path, nfof_lim , color, name_tag):
        self.name = name
        self.f_sat = []
        self.n_tot = []
        self.mean_mass = []
        self.path = path
        
        self.color = {'all': color , 'cen': color  , 'sat': color  }
             
        self.name_tag =  name_tag      
        self.mass_arr= mass_arr
        self.nfof_lim = nfof_lim
        self.thickness = 1

    #these plotting options are different for the different sub objects    

    def redefine_colors(self,dic):
        self.color = dic

    def redefine_linestyles(self,dic):
        self.linestyle = dic
        
    def redefine_thickness(self,dic):
        self.thickness = dic

    def redefine_markers(self, dic):
        self.marker = dic    

class ESD_sim(ESD):
    def __init__(self,  name , mass_arr, path , nfof_lim, color, name_tag, m_host_limit, extra_tag):
            ESD.__init__(self,  name , mass_arr, path , nfof_lim, color, name_tag)
            self.m_host_limit = m_host_limit
            self.extra_tag =  extra_tag
            marker = ''
            self.marker = {'all': marker  , 'cen':  marker  , 'sat': marker  }   
            linestyle = '-'
            self.linestyle = {'all': linestyle  , 'cen': linestyle   , 'sat': linestyle  }   
            self.HUBBLEPARAM = 0.677700 # this should be read from the file
            self.sub_type = ['all' , 'sat' ,'cen']
            filename_sim = self.get_filename(0, 'TOT')
            ESD_bins = 150
            self.ESD_bins = ESD_bins
            
            
            dt_ESD_profiles = np.dtype([('mass_bin', str , 2) , \
                                    ('rad' , float , ESD_bins  ), \
                                    ('all' , float , ESD_bins  ), \
                                    ('sat' , float , ESD_bins  ), \
                                    ('cen' , float , ESD_bins  ) ])
            
        
            self.ESD_profiles =np.zeros(len(self.mass_arr)-1 , dtype=dt_ESD_profiles)
            #84th percentile - mean for plotting pourpose
            self.ESD_84th =np.zeros(len(self.mass_arr)-1 , dtype=dt_ESD_profiles)
            #mean - 16th percentile for plotting pourpose
            self.ESD_16th =np.zeros(len(self.mass_arr)-1 , dtype=dt_ESD_profiles)


    def ESD_in_bins(self):
        """ computes the ESD profile by averaging the ESD profiles of the galaxies in the bin
        """
        n_ESD_bins = self.ESD_bins
        dt_ESD_template = np.dtype([('mean_snap_esd_profile_y'  , float , n_ESD_bins), \
                                    ('stdev_snap_esd_profile_y'  , float , n_ESD_bins), \
                                    ('boot_16th_snap_esd_profile_y'  , float , n_ESD_bins), \
                                    ('boot_84th_snap_esd_profile_y'  , float , n_ESD_bins), \
                                    ('snap_r_plus_half' , float , n_ESD_bins )   ])
        ESD_prof_in_bin  = { 'all' : np.zeros(1, dtype=dt_ESD_template) , \
                             'cen' : np.zeros(1, dtype=dt_ESD_template) , \
                             'sat' : np.zeros(1, dtype=dt_ESD_template)}
        self.ESD_prof_in_bin = []
        for ibin in np.arange(len(self.mass_arr)-1) :
            #print ibin , ', bin ' , self.mass_arr[ibin] , self.mass_arr[ibin+1]
                
            for sub_type_tag   in  self.sub_type:

                if sub_type_tag == 'cen':
                    index_select , = where((self.EAGLE_gal_cat['selected'][:] == ibin+1) & \
                                           (self.EAGLE_gal_cat['main_subhalo'][:] == 1 ))
                elif sub_type_tag == 'sat':
                    index_select , = where((self.EAGLE_gal_cat['selected'][:] == ibin+1) & \
                                           (self.EAGLE_gal_cat['main_subhalo'][:] == 0 ))
                else :
                    index_select , = where((self.EAGLE_gal_cat['selected'][:] == ibin+1))

                #print 'selected subhaloes' , len(index_select)
                if  len(index_select)  ==0 :
                    print 'empty bin'
                    continue

                
                ESD_prof_in_bin[sub_type_tag]['snap_r_plus_half'][0,:] = self.EAGLE_gal_cat['snap_r_plus_half'][index_select[0]][:]
                #print np.shape(ESD_prof_in_bin[sub_type_tag]['mean_snap_esd_profile_y'][0,:])

                mean_esd =np.zeros(n_ESD_bins, float)

                for i in np.arange(n_ESD_bins):
                    value=np.sum(self.EAGLE_gal_cat['snap_esd_profile_y'][index_select,i])                   
                    
                    mean_esd[i] = value/len(index_select)

                        

                ESD_prof_in_bin[sub_type_tag]['mean_snap_esd_profile_y'][0,:] = mean_esd/100.

                n_realiz =1000
                seed = (np.random.rand(n_realiz)*10.0**5).astype(int64)
                mean_esd =np.zeros((n_realiz,n_ESD_bins), float)
                for b in np.arange(n_realiz) :
                    np.random.seed(seed=seed[b])                  
                    random_indexes =  np.floor(np.random.rand(len(index_select))*len(index_select)).astype(int)
                    #print random_indexes 
                    index_select_rnd = index_select[random_indexes]
                    #print index_select_rnd                     
                    


                    for i in np.arange(n_ESD_bins):
                        value=np.sum(self.EAGLE_gal_cat['snap_esd_profile_y'][index_select_rnd,i])                   
                        
                        mean_esd[b,i] = value/len(index_select)

                        

                
                for i in np.arange(n_ESD_bins):
                    ESD_prof_in_bin[sub_type_tag]['boot_16th_snap_esd_profile_y'][0,i] = np.percentile(mean_esd[:,i],16)/100.
                    
                    ESD_prof_in_bin[sub_type_tag]['boot_84th_snap_esd_profile_y'][0,i] = np.percentile(mean_esd[:,i],84)/100.

                
                #print vec
                
            temp = cp.deepcopy(ESD_prof_in_bin)
            self.ESD_prof_in_bin.append(temp)


                
            

    def stat_on_bins(self):
        """ it computes the statistics in the bins and save them in a record file"""
        # ESD_object.vec_stat[0]['cen']['histo_y_subhalo_mass']
        #consider to save histograms for imp quantities in bin, stellar mass, halo mass
        # luminosity color
        # save median masses for halo and stars

        
        #find the list of subhaloes
        histo_bins=15
        self.stat_histo_bins = histo_bins

        dt_histo_template = np.dtype([('mean_subhalo_mass', float ) , \
                                      ('median_subhalo_mass', float ) , \
                                      ('16th_p_subhalo_mass', float ) , \
                                      ('84th_p_subhalo_mass', float ) , \
                                      ('histo_x_subhalo_mass' , float , histo_bins  ), \
                                      ('histo_y_subhalo_mass' , float , histo_bins  ), \
                                      ('mean_subhalo_half_mass_rad_dm', float ) , \
                                      ('median_subhalo_half_mass_rad_dm', float ) , \
                                      ('16th_p_subhalo_half_mass_rad_dm', float ) , \
                                      ('84th_p_subhalo_half_mass_rad_dm', float ) , \
                                      ('histo_x_subhalo_half_mass_rad_dm' , float , histo_bins  ), \
                                      ('histo_y_subhalo_half_mass_rad_dm' , float , histo_bins  ), \
                                      ('mean_main_halo_m200_crit', float ) , \
                                      ('median_main_halo_m200_crit', float ) , \
                                      ('16th_p_main_halo_m200_crit', float ) , \
                                      ('84th_p_main_halo_m200_crit', float ) , \
                                      ('histo_x_main_halo_m200_crit' , float , histo_bins  ), \
                                      ('histo_y_main_halo_m200_crit' , float , histo_bins  ), \
                                      ('mean_sub_halo_rad_dist_main', float ) , \
                                      ('median_sub_halo_rad_dist_main', float ) , \
                                      ('16th_p_sub_halo_rad_dist_main', float ) , \
                                      ('84th_p_sub_halo_rad_dist_main', float ) , \
                                      ('histo_x_sub_halo_rad_dist_main' , float , histo_bins  ), \
                                      ('histo_y_sub_halo_rad_dist_main' , float , histo_bins  ), \
                                      ('mean_subhalo_mass_in_star', float ) , \
                                      ('median_subhalo_mass_in_star', float ) , \
                                      ('16th_p_subhalo_mass_in_star', float ) , \
                                      ('84th_p_subhalo_mass_in_star', float ) , \
                                      ('histo_x_subhalo_mass_in_star' , float , histo_bins  ), \
                                      ('histo_y_subhalo_mass_in_star' , float , histo_bins  ), \
                                      ('mean_subhalo_mass_in_star_30kpc', float ) , \
                                      ('median_subhalo_mass_in_star_30kpc', float ) , \
                                      ('16th_p_subhalo_mass_in_star_30kpc', float ) , \
                                      ('84th_p_subhalo_mass_in_star_30kpc', float ) , \
                                      ('histo_x_subhalo_mass_in_star_30kpc' , float , histo_bins  ), \
                                      ('histo_y_subhalo_mass_in_star_30kpc' , float , histo_bins  ), \
                                      ('mean_Msub_Mhost200_ratio', float ) , \
                                      ('median_Msub_Mhost200_ratio', float ) , \
                                      ('16th_p_Msub_Mhost200_ratio', float ) , \
                                      ('84th_p_Msub_Mhost200_ratio', float ) , \
                                      ('histo_x_Msub_Mhost200_ratio' , float , histo_bins  ), \
                                      ('histo_y_Msub_Mhost200_ratio' , float , histo_bins  ), \
                                     
                                      ('n_sub', int )  \
        ])

        
                                     
        dt_histo_mass = np.dtype([('mean_value', float ) , \
                                  ('histo_mass_bin' , float , histo_bins  ), \
                                  ('subhalo_mass_all' , float , histo_bins ), \
                                  ('sat_host_mass' , float ,  histo_bins  ), \
                                  ('cen_m200_crit' , float , histo_bins   ) ])
        
        self.histo_halo_mass =np.zeros(len(self.mass_arr)-1 , dtype=dt_histo_mass)
        self.avail_stat = ['subhalo_mass' , \
                           'subhalo_mass_in_star' ,\
                           'subhalo_mass_in_star_30kpc', \
                           'main_halo_m200_crit' , \
                           'subhalo_half_mass_rad_dm', \
                           'sub_halo_rad_dist_main', \
                           'Msub_Mhost200_ratio']

        stat_subhalo_mass  = { 'all' : np.zeros(1, dtype=dt_histo_template) , \
                               'cen' : np.zeros(1, dtype=dt_histo_template) , \
                               'sat' : np.zeros(1, dtype=dt_histo_template)}
        self.vec_stat = []
        for ibin in np.arange(len(self.mass_arr)-1) :
            print ibin , ', bin ' , self.mass_arr[ibin] , self.mass_arr[ibin+1]
                
            for sub_type_tag   in  ['all' , 'sat' ,'cen']:

                if sub_type_tag == 'cen':
                    index_select , = where((self.EAGLE_gal_cat['selected'][:] == ibin+1) & \
                                           (self.EAGLE_gal_cat['main_subhalo'][:] == 1 ))
                elif sub_type_tag == 'sat':
                    index_select , = where((self.EAGLE_gal_cat['selected'][:] == ibin+1) & \
                                           (self.EAGLE_gal_cat['main_subhalo'][:] == 0 ))
                else :
                    index_select , = where((self.EAGLE_gal_cat['selected'][:] == ibin+1))

                stat_subhalo_mass[sub_type_tag]['n_sub'] = len(index_select)
                if  stat_subhalo_mass[sub_type_tag]['n_sub'] ==0 :
                    print 'empty bin'
                    continue
                for mass_type_tag in self.avail_stat:
                    
                    
                    if mass_type_tag == 'Msub_Mhost200_ratio' :
                       vec_mass = np.log10(10.0**self.EAGLE_gal_cat['subhalo_mass'][index_select] / 10.0**self.EAGLE_gal_cat['main_halo_m200_crit'][index_select])
                    elif (mass_type_tag == 'subhalo_half_mass_rad_dm') | (mass_type_tag == 'sub_halo_rad_dist_main'):
                       vec_mass = np.log10(self.EAGLE_gal_cat[mass_type_tag][index_select])
                       
                    else:
                       vec_mass = self.EAGLE_gal_cat[mass_type_tag][index_select]
                        
                    stat_subhalo_mass[sub_type_tag]['mean_'+mass_type_tag] = np.log10(mean(10.0**vec_mass))
                    stat_subhalo_mass[sub_type_tag]['median_'+mass_type_tag] = np.log10(median(10.0**vec_mass))
                    stat_subhalo_mass[sub_type_tag]['16th_p_'+mass_type_tag] = np.log10(np.percentile(10.0**vec_mass,16))
                    stat_subhalo_mass[sub_type_tag]['84th_p_'+mass_type_tag] = np.log10(np.percentile(10.0**vec_mass,84))
                    #print sub_type_tag , ' ' , mass_type_tag , ' mean ' , stat_subhalo_mass[sub_type_tag]['mean_'+mass_type_tag], ' min max ' , min(vec_mass) , max(vec_mass), 'median' , stat_subhalo_mass[sub_type_tag]['median_'+mass_type_tag] , '16 84 ' ,  stat_subhalo_mass[sub_type_tag]['16th_p_'+mass_type_tag]  ,  stat_subhalo_mass[sub_type_tag]['84th_p_'+mass_type_tag] 
                    if (min(vec_mass) < 1) &  (mass_type_tag == 'subhalo_mass') : print '*****************************************'
                    bin_w =  (self.mass_arr[ibin+1]- self.mass_arr[ibin])

                    if (mass_type_tag == 'Msub_Mhost200_ratio') | (mass_type_tag == 'subhalo_half_mass_rad_dm') | (mass_type_tag == 'sub_halo_rad_dist_main'):
                        min_h = min(vec_mass)    
        
                    else:
                        if sum(vec_mass) ==0 :
                            min_h =0
                        else:
                            min_h = min(vec_mass[vec_mass>0])
                            
                            
                    max_h = max(vec_mass)
                    #print min_h , max_h
                    hist, bin_edges = np.histogram(vec_mass, bins=size(stat_subhalo_mass[sub_type_tag]['histo_y_'+mass_type_tag]), range=(min_h , max_h),  density=True)

                    #print size(stat_subhalo_mass[sub_type_tag]['histo_x_'+mass_type_tag][0,:])
                    stat_subhalo_mass[sub_type_tag]['histo_x_'+mass_type_tag][0,:] = [(bin_edges[i] + bin_edges[i+1])/2. for i in range(len(bin_edges)-1)]
                    stat_subhalo_mass[sub_type_tag]['histo_y_'+mass_type_tag][0,:] = hist[:]                  
                                                                                           
           
            
            temp = cp.deepcopy(stat_subhalo_mass)
            #print stat_subhalo_mass['all']['histo_x_subhalo_mass']
            #print stat_subhalo_mass['all']['histo_y_subhalo_mass']
            self.vec_stat.append(temp)

            
        #print self.vec_stat['all']['mean_subhalo_mass'][2,:]

        #print self.vec_stat[0]['all']['histo_x_subhalo_mass']
        #print self.vec_stat[1]['all']['histo_x_subhalo_mass']        
        #print self.vec_stat[2]['all']['histo_x_subhalo_mass']

       

            
    def read_EAGLE_gal_cat(self):

        """ read the ESD profiles as computed by the IDL code"""


        #get the *.fits file name for the given simulation
        # apply the same halo selection in stellar mass limit
        path = self.path.replace("REFERENCE_PLANCK_L0100N1504_ESD_stat_Z_MsunM_host_limit_","")

        
        
        sim_cat = path + 'all_sub_prop_EAGLE_sub8p0.fits'

        print 'open ' + sim_cat
        #it is not opening the header
        #hdulist = pyfits.open(ShearMergedCatalogueGroups)[1].data
        hdulist = pyfits.open(sim_cat)[1].data
        mass = hdulist.field("subhalo_mass")
        n_el_EAGLE_cat = len(mass)
        n_bins=self.ESD_bins


        #create the catalogue structure with recarrays
        dt_EAGLE_gal_cat = np.dtype([('subhalo_id', long), \
                                     ('subhalo_mass_in_star', float), \
                                     ('main_halo_m200_crit', float), \
                                     ('subhalo_mass', float), \
                                     ('subhalo_mass_in_star_30kpc', float), \
                                     ('subhalo_half_mass_rad_dm', float), \
                                     ('sub_halo_rad_dist_main', float), \
                                     ('ih_parent_fof' , long), \
                                     ('main_subhalo' , int), \
                                     ('selected' , int ), \
                                     ('galaxyID' , long), \
                                     ('snap_r_plus_half' , float , n_bins ), \
                                     ('snap_esd_profile_y'  , float , n_bins )   ])

        self.EAGLE_gal_cat =np.zeros(n_el_EAGLE_cat , dtype=dt_EAGLE_gal_cat)

        self.EAGLE_gal_cat['subhalo_mass'][:] = hdulist.field("subhalo_mass")        
        self.EAGLE_gal_cat['subhalo_id'][:] = hdulist.field("subhalo_id")
        self.EAGLE_gal_cat['subhalo_mass_in_star'][:] = hdulist.field("subhalo_mass_in_star")
        self.EAGLE_gal_cat['subhalo_mass_in_star_30kpc'][:] = hdulist.field("subhalo_mass_in_star_30kpc")
        self.EAGLE_gal_cat['ih_parent_fof'][:] = hdulist.field("ih_parent_fof")
        self.EAGLE_gal_cat['main_subhalo'][:] = hdulist.field("main_subhalo")
        self.EAGLE_gal_cat['main_halo_m200_crit'][:] = hdulist.field("main_halo_m200_crit")
        self.EAGLE_gal_cat['subhalo_half_mass_rad_dm'][:] = hdulist.field("subhalo_half_mass_rad_dm")
        self.EAGLE_gal_cat['sub_halo_rad_dist_main'][:] = hdulist.field("sub_halo_rad_dist_main")
        self.EAGLE_gal_cat['snap_esd_profile_y'][:] = hdulist.field("snap_esd_profile_y")
        self.EAGLE_gal_cat['snap_r_plus_half'][:] = hdulist.field("snap_r_plus_half") 
        #self.EAGLE_gal_cat['galaxyID'][:] = hdulist.field("galaxyID")
        #self.EAGLE_gal_cat['color'][:] = hdulist.field("color")
        #uncomment it when using the matched cat



        #check if the selection is done in stellar mass or stellar mass in 30kp
        if self.extra_tag == 'Mstar30kpc' :            
            vector_sub_mass =  self.EAGLE_gal_cat['subhalo_mass_in_star_30kpc'][:]
        else:
            vector_sub_mass =  self.EAGLE_gal_cat['subhalo_mass_in_star'][:]

            
        if min(vector_sub_mass) > min(self.m_host_limit) :
            raise NameError('catalogue does not contain galaxies with masses lower than the M_lim')


            
        for ibin in np.arange(len(self.mass_arr)-1) :
            

            nFOF_limit = self.nfof_lim

            HUBBLEPARAM =self.HUBBLEPARAM
            mass_star_limit_fof = self.m_host_limit[ibin]  + log10(HUBBLEPARAM) #!!!!!!!! check

            print ibin , ', bin ' , self.mass_arr[ibin] , self.mass_arr[ibin+1] 
            print 'm host limit ' , self.m_host_limit[ibin] , ' nfof limit ' , nFOF_limit
            
            # list of sub over the mass limit
            index_sub_over_mstar_limit ,  = np.where(vector_sub_mass > mass_star_limit_fof)
            ind_bool_sub_over_mstar_limit = vector_sub_mass > mass_star_limit_fof
            n_el_over_mstar_limit = sum(ind_bool_sub_over_mstar_limit)
            n_index_sub_over_mstar_limit = len(index_sub_over_mstar_limit)
            #print 'n el above the mstar limit' , n_el_over_mstar_limit ,  len(index_sub_over_mstar_limit)

            
            
            # highest FOF number of the sub over the mass limit
            max_ih_parent_fof  = max(self.EAGLE_gal_cat['ih_parent_fof'][index_sub_over_mstar_limit])

            #print 'max ind parent FOF' , max_ih_parent_fof

            #number of subhalos for every fof group that are above the mass limit
            N_sub_over_limit_parent_fof  =  np.zeros(max_ih_parent_fof+1, int)
            
            for  i  in  arange(n_index_sub_over_mstar_limit)  :
                N_sub_over_limit_parent_fof[self.EAGLE_gal_cat['ih_parent_fof'][index_sub_over_mstar_limit[i]]] +=1
            print 'total subhaloes above the mass limitn in FOF groups' , sum(N_sub_over_limit_parent_fof)

            ih_fof_over_limit, = np.where(N_sub_over_limit_parent_fof >= nFOF_limit)
            n_fof_over_limit = len(ih_fof_over_limit)

            print n_fof_over_limit  , ' FOF froups with more than ' ,nFOF_limit , 'members above the mass limit ' , mass_star_limit_fof
      # ;list of indexes of FOF over the limit
                
            #  compute the number of sub for every FOF

            
            bool_of_subhaloes_in_rich_goups = np.zeros(n_el_EAGLE_cat, long) 
            for  i  in  arange(n_fof_over_limit) :
                index,  = np.where(self.EAGLE_gal_cat['ih_parent_fof'][:] ==  ih_fof_over_limit[i])
                bool_of_subhaloes_in_rich_goups[index] += 1
            

            #print max(bool_of_subhaloes_in_rich_goups)
            min_mass_bin = self.mass_arr[ibin] + log10(HUBBLEPARAM)
            max_mass_bin = self.mass_arr[ibin+1]+ log10(HUBBLEPARAM)

            index_selected_rich_sub_in_mass_bin , = where( (vector_sub_mass > min_mass_bin) & (vector_sub_mass <=  max_mass_bin) & ( bool_of_subhaloes_in_rich_goups > 0) )

            self.EAGLE_gal_cat['selected'][index_selected_rich_sub_in_mass_bin] = ibin +1

            #print  'number of selected sub ' , len(index_selected_rich_sub_in_mass_bin)

        for ibin in np.arange(len(self.mass_arr)-1) :
                index, = where( self.EAGLE_gal_cat['selected'][:] == ibin+1)
                print 'bin ', ibin+1 , 'contain ' , len(index) , ' galaxies'
            
            
    def compute_chi_squared(self, ESD_gama, sub_type_tag, path):
            """computes the chi square between model and data"""
            chi_sq = 0
            min_value = []
            sim = []
            data = []
            dof=0.
            #print ESD_gama.ESD_profiles[0]['rad'][:]

            filename =path+'_chi_squared_'+ self.name +'_'+ sub_type_tag+'.txt'
            with open(filename, "w") as text_file:
                text_file.write(" bin_left ,"+\
                                 " bin_right ,"+\
                                 " rad ,"+\
                                 " chi_sq_term \n")
                for ibin  in np.arange(1,len( self.mass_arr)-1):


                        for i  in np.arange(len( ESD_gama.ESD_profiles[ibin]['rad'][:])):

                                # index at which the distance between the sim and gama is at min

                                min_index = np.argmin(abs(self.ESD_profiles[ibin]['rad'] - ESD_gama.ESD_profiles[ibin]['rad'][i]))
                                #print  min_index
                                sim.append(self.ESD_profiles[ibin][sub_type_tag][min_index] )
                                data.append(ESD_gama.ESD_profiles[ibin][sub_type_tag][i] )
                                # chi_sq
                                # (esd_gama - esd_sim)^2 / (sigma_gama^2 + sigma_sim^2)
                                chi_sq_term =   (self.ESD_profiles[ibin][sub_type_tag][min_index] - \
                                             ESD_gama.ESD_profiles[ibin][sub_type_tag][i])**2. / \
                                (
                                ( (ESD_gama.ESD_84th[ibin][sub_type_tag][i]-
                                   ESD_gama.ESD_16th[ibin][sub_type_tag][i])/2 )**2 +

                                ( (self.ESD_84th[ibin][sub_type_tag][min_index]  - \
                                   self.ESD_16th[ibin][sub_type_tag][min_index])/2 )**2)
                                chi_sq += chi_sq_term
                                text_file.write((" {:03.2f} ,"+\
                                                 " {:03.2f} ,"+\
                                                 " {:03.4f} ,"+\
                                                 " {:03.4f} \n")\
                                                .format(self.mass_arr[ibin],\
                                                        self.mass_arr[ibin+1],\
                                                        ESD_gama.ESD_profiles[ibin]['rad'][i],\
                                                        chi_sq_term))

                                dof +=1.


                return chi_sq / (dof-1) ,  (dof-1)
            
    def get_ESD_profile(self):
        """ collects the esd profiles computed by the ESD_in_bins function and puts them in the proper
        container (same as the KIDS data containers)"""
        sub_type = ['all' , 'sat' ,'cen']
        sub_type_sim = ['TOT', 'SAT' , 'CEN']

        # if the data is computed use it otherwise it is looking for the IDL saved files
        if hasattr(self, 'ESD_prof_in_bin') :
            for type_index, sub_type_tag in enumerate(sub_type):
                    for ibin in np.arange( len(self.mass_arr)-1) :
                        
                            self.ESD_profiles[ibin][sub_type_tag][:] = self.ESD_prof_in_bin[ibin][sub_type_tag]['mean_snap_esd_profile_y'][0,:]
                            self.ESD_profiles[ibin]['rad'][:] = self.ESD_prof_in_bin[ibin][sub_type_tag]['snap_r_plus_half'][0,:]
                            self.ESD_16th[ibin][sub_type_tag][:] = (self.ESD_prof_in_bin[ibin][sub_type_tag]['boot_16th_snap_esd_profile_y'][0,:])

                            # !!! this is for visualization only ESD_84th_error is the portion of the error bar that goues from the mean to the top error
                            self.ESD_84th[ibin][sub_type_tag][:]= (self.ESD_prof_in_bin[ibin][sub_type_tag]['boot_84th_snap_esd_profile_y'][0,:] )
                            self.ESD_84th[ibin]['rad'][:] = self.ESD_profiles[ibin]['rad'][:]
                            self.ESD_16th[ibin]['rad'][:] = self.ESD_profiles[ibin]['rad'][:]
                        

        else :
            print 'compute the ESD first'

            # print self.mass_arr
            # for type_index, sub_type_tag in enumerate(sub_type):
            #         for ibin in np.arange( len(self.mass_arr)-1) :
            #                 filename_sim = self.get_filename(ibin, sub_type_sim[type_index])


            #                 radius, mean_ESD , median_ESD , ESD_down_std, ESD_up_std ,  ESD_error , ESD_16th_down_error , ESD_84th_up_error , ESD_down_error , ESD_up_error = loadtxt(filename_sim,unpack=True, comments='#')

            #                 self.ESD_profiles[ibin][sub_type_tag][:] = pow(10.0 , mean_ESD[:])
            #                 self.ESD_profiles[ibin]['rad'][:] = 10.0**radius[:]
            #                 self.ESD_16th[ibin][sub_type_tag][:] = -10.0**ESD_down_error[:] + 10.0**mean_ESD[:]

            #                 self.ESD_84th[ibin][sub_type_tag][:]= 10.0**ESD_up_error[:] - 10.0**mean_ESD[:]
            #                 self.ESD_16th_error[ibin]['rad'][:] = 10.0**radius[:]
            #                 self.ESD_84th_error[ibin]['rad'][:] = 10.0**radius[:]







    def get_filename(self, ibin, type_tag):                        
            min_mass_tag = '%.1f' %  self.mass_arr[ibin]
            max_mass_tag = '%.1f' %  self.mass_arr[ibin+1]
                            
            sim_mass_limit = '%.2f' %  self.m_host_limit[ibin]
            #print sim_mass_limit  + ' ' +min_mass_tag+ ' ' +max_mass_tag
            min_mass_tag_with_p = min_mass_tag.replace(".", "p")
                        
            max_mass_tag_with_p =  max_mass_tag.replace(".", "p")
            filename_sim = self.path + sim_mass_limit + 'Nfof_gt'+  str(self.nfof_lim-1)+ '_min_mass_' + min_mass_tag_with_p+ '__max_mass_'+ max_mass_tag_with_p+'__max_rad_2p0boot1000' + self.extra_tag +type_tag + '.dat'
            #print  filename_sim
            return filename_sim
                            
                           
    def get_fsat(self):
            """get the f_sat info from the data"""
            f_sat = []
            mean_mass = []

            if hasattr(self, 'vec_stat') :
                for ibin in np.arange( len(self.mass_arr)-1) :

                    value = float(self.vec_stat[ibin]['sat']['n_sub'][0])/ float(self.vec_stat[ibin]['all']['n_sub'][0])
                    f_sat.append( value   )       
                    mean_mass.append( (self.mass_arr[ibin] + self.mass_arr[ibin+1] )/2.0 )        
            else:



            
                for ibin in np.arange( len(self.mass_arr)-1) :
                        filename_sim = self.get_filename(ibin, 'TOT')

                        value = []

                        #the with statement opne and close the file automatically 
                        with open(filename_sim, 'r') as f:
                                first_line = f.readline()
                                second_line = f.readline()

                        #print second_line

                        #for every space separated word it tries to append the float conversion of it to a value (that is a list)
                        # if the try is failing it will pass to the next iteration
                        for t in second_line.split():
                                try:
                                        value.append(float(t))
                                except ValueError:
                                        pass

                        # here since there is only one number I assume that the previus code will return one number only and not a list
                        #print value

                        f_sat.append( 1 - value[0]   )       
                        mean_mass.append( (self.mass_arr[ibin] + self.mass_arr[ibin+1] )/2.0 )


            self.f_sat = f_sat
            self.mean_mass = mean_mass





            


class ESD_gama(ESD):
    def __init__(self, name , mass_arr , path , nfof_lim, color, name_tag, extra_tag , cat_blind , m_bias):
            ESD.__init__(self,name , mass_arr,  path, nfof_lim , color, name_tag)
            marker = 'o'
            linestyle = ''
            #default linestyle
            self.linestyle = {'all': linestyle  , 'cen': linestyle   , 'sat': linestyle  }
            #default markers
            self.marker = {'all': marker  , 'cen': marker    , 'sat': marker   }
            #map between file name and type of galaxies
            self.dict_type={ 'sat' : 'RankBCG2-inf' , 'cen' : 'RankBCG1' , 'all' : 'RankBCG1-inf'}
            #multiplicative bias
            self.m_bias = m_bias
            self.m_bias_corr  = 1. / (1. + m_bias)
            self.extra_tag =  extra_tag
            # blind catalogue's name
            self.cat_blind = cat_blind

            #compute the number of bins in the file by counting the number of lines
            ESD_bins = sum(1 for line in open(self.get_filename(0, 'all')))-1
            self.ESD_bins = ESD_bins
            print self.ESD_bins

            #datatype for the record array
            dt_ESD_profiles = np.dtype([('mass_bin', str , 2) , \
                                        ('rad' , float , ESD_bins  ), \
                                        ('all' , float , ESD_bins  ), \
                                        ('sat' , float , ESD_bins  ), \
                                        ('cen' , float , ESD_bins  ) ])
            
            #stored record array with the data and the errors
            self.ESD_profiles =np.zeros(len(self.mass_arr)-1 , dtype=dt_ESD_profiles)
            self.ESD_16th =np.zeros(len(self.mass_arr)-1 , dtype=dt_ESD_profiles)
            self.ESD_84th =np.zeros(len(self.mass_arr)-1 , dtype=dt_ESD_profiles)


    def get_filename(self, ibin, type_tag):
        """ function that computes the filename for the KIDS data given galaxy type
        mass bin
        Nfof
        and the extra_tag (the number of galaxies of the given selection
        !!! watch out that there is no relation now between the bin number and the mass bin
        !!! this has to be taken care by the user
        """
        dict={ 'sat' : 'RankBCG2-inf' , 'cen' : 'RankBCG1' , 'all' : 'RankBCG1-inf'}

        filename = self.path+'shearcovariance_logmstarbin' + str(ibin+1) + 'of' + str(len(self.mass_arr)-1)+'_'+ self.extra_tag+'_'  + 'Nfof' + str(self.nfof_lim)+'-inf'+'_'+ self.dict_type[type_tag]  +  '_'+ self.cat_blind + '.txt'
        return filename

    def get_filename_lensIDs(self, ibin, type_tag):
        """ function that returns the filename where the lenses id are stored
        this is used to compute the fsat from the data"""
        dict={ 'sat' : 'RankBCG2-inf' , 'cen' : 'RankBCG1' , 'all' : 'RankBCG1-inf'}

        filename = self.path+'shearcovariance_logmstarbin' + str(ibin+1) + 'of' + str(len(self.mass_arr)-1)+'_'+ self.extra_tag+'_'  + 'Nfof' + str(self.nfof_lim)+'-inf'+'_'+ self.dict_type[type_tag]  +  '_'+ 'lensIDs.txt'
        return filename


        
    def get_ESD_profile(self):
            """ load the ESD profile given by the KIDS pipeline 
            same warning for the bin number as for the get_filename func
            """
            
            sub_type = ['all' , 'sat' ,'cen']
            
            for sub_type_tag in sub_type:
                    for massbin in np.arange(len(self.mass_arr)-1):
                            
                            min_mass_tag = '%.1f' %  self.mass_arr[massbin]
                            max_mass_tag = '%.1f' %  self.mass_arr[massbin+1]
                            
                            bin_tag = '%i' % (massbin + 2) #data starts from bin 1 -> 9.7,10

                            filename = self.get_filename(massbin, sub_type_tag)
                            
                            radius , ESD , aaaa , error_ESD , a , b , c , d = loadtxt(filename,unpack=True,comments='#')
                            
                            self.ESD_profiles[massbin]['mass_bin'] = [min_mass_tag ,  \
                                                                      max_mass_tag]
                            self.ESD_profiles[massbin]['rad'][:]= radius[:]/1000.
                            self.ESD_profiles[massbin][sub_type_tag][:]= ESD[:] * self.m_bias_corr
                            self.ESD_16th[massbin]['rad'][:]= radius[:]/1000.
                            self.ESD_16th[massbin][sub_type_tag][:]= self.ESD_profiles[massbin][sub_type_tag][:] - error_ESD[:]* self.m_bias_corr
                            self.ESD_84th[massbin]['rad'][:]= radius/1000.

                            
                            self.ESD_84th[massbin][sub_type_tag][:]= self.ESD_profiles[massbin][sub_type_tag][:] + error_ESD[:]* self.m_bias_corr

    def get_fsat(self):
            "get the f_sat info from the data"
            f_sat = []
            mean_mass = []
            for ibin in np.arange(len(self.mass_arr)-1) :
                    

                    filename_sat = self.get_filename_lensIDs(ibin, 'sat')
                    filename_all = self.get_filename_lensIDs(ibin, 'all')
                    
                    n_sat = sum(1 for line in open(self.get_filename_lensIDs(ibin, 'sat')))-1
                    n_tot = sum(1 for line in open(self.get_filename_lensIDs(ibin, 'all')))-1
                    #print filename
                    #print filename_sat
                    #print filename_all
                    
                    x = float(n_sat) / float(n_tot)
                    f_sat.append(x)
                    mean_mass.append( (self.mass_arr[ibin] + self.mass_arr[ibin+1])/2.0 ) 
                    #print self.mass_arr[ibin] , self.mass_arr[ibin+1] , n_sat , n_tot , x
                    
            self.f_sat = f_sat
            self.mean_mass = mean_mass







def plot_ESD_allbins(ESD_data, savename):
    """ plots the esd in every bin for the objects in ESD_data """
    
                        
    scale_color = 40
    sub_type = ['all' , 'sat' ,'cen']

    sub_type_cm = [plt.cm.gray , plt.cm.winter ,plt.cm.autumn]
    #----------- plotting the profiles in diff bins

    index_50kpc =np.argmin(abs(ESD_data.ESD_profiles['rad'][0,:] - 0.05))
    #print index_50kpc , ESD_data.ESD_profiles['rad'][0,index_50kpc]
    index_500kpc =np.argmin(abs(ESD_data.ESD_profiles['rad'][0,:] - 0.5))
    #print index_500kpc , ESD_data.ESD_profiles['rad'][0,index_500kpc]
    
    if isinstance(ESD_data, ESD_sim):

        for sub_type_tag , color_cm  in zip(sub_type , sub_type_cm):

            #print np.log10([10.0**((ESD_data.mass_arr[ibin]+ESD_data.mass_arr[ibin+1])/2.0)  for ibin in  np.arange(len(ESD_data.mass_arr)-1)])
            #print np.log10([(ESD_data.ESD_profiles[sub_type_tag][ibin,index_50kpc]/ESD_data.ESD_profiles[sub_type_tag][0,index_50kpc]) for ibin in  np.arange(len(ESD_data.mass_arr)-1)])
            #print ESD_data.ESD_profiles['rad'][0,:]


            plt.figure(figsize=(7,6))
            plt.loglog([10.0**((ESD_data.mass_arr[ibin]+ESD_data.mass_arr[ibin+1])/2.0) for ibin in  np.arange(len(ESD_data.mass_arr)-1)] ,\
                       [(ESD_data.ESD_profiles[sub_type_tag][ibin,index_50kpc]/ESD_data.ESD_profiles[sub_type_tag][0,index_50kpc]) for ibin in  np.arange(len(ESD_data.mass_arr)-1)] ,\
                       color =color_cm(0), \
                       label= '$\Delta \Sigma^{'+ sub_type_tag +'}(r_p = 0.05Mpc) $ ',\
                       marker=ESD_data.marker[sub_type_tag] ,\
                       ls =ESD_data.linestyle[sub_type_tag])


            plt.loglog([10.0**((ESD_data.mass_arr[ibin]+ESD_data.mass_arr[ibin+1])/2.0) for ibin in  np.arange(len(ESD_data.mass_arr)-1)] ,\
                       [10.0**(ESD_data.vec_stat[ibin][sub_type_tag]['mean_subhalo_mass']-ESD_data.vec_stat[0][sub_type_tag]['mean_subhalo_mass']) for ibin in  np.arange(len(ESD_data.mass_arr)-1)] ,\
                       color ='k', \
                       label= '$M_{sub}^{'+ sub_type_tag +'}$ ',\
                       marker=ESD_data.marker[sub_type_tag] ,\
                       ls =ESD_data.linestyle[sub_type_tag])


            if sub_type_tag == 'sat' :
                plt.loglog([10.0**((ESD_data.mass_arr[ibin]+ESD_data.mass_arr[ibin+1])/2.0) for ibin in  np.arange(len(ESD_data.mass_arr)-1)] ,\
                           [(ESD_data.ESD_profiles[sub_type_tag][ibin,index_500kpc]/ESD_data.ESD_profiles[sub_type_tag][0,index_500kpc]) for ibin in  np.arange(len(ESD_data.mass_arr)-1)] ,\
                           color =color_cm(0), \
                           label= '$\Delta \Sigma^{'+ sub_type_tag +'}(r_p = 0.5Mpc) $ ',\
                           marker=ESD_data.marker[sub_type_tag] ,\
                            ls ='--')


                plt.loglog([10.0**((ESD_data.mass_arr[ibin]+ESD_data.mass_arr[ibin+1])/2.0) for ibin in  np.arange(len(ESD_data.mass_arr)-1)] ,\
                           [10.0**(ESD_data.vec_stat[ibin][sub_type_tag]['mean_main_halo_m200_crit']-ESD_data.vec_stat[0][sub_type_tag]['mean_main_halo_m200_crit']) for ibin in  np.arange(len(ESD_data.mass_arr)-1)] ,\
                           color ='k', \
                           label= '$M_{sub}^{'+ sub_type_tag +'}$ ',\
                           marker=ESD_data.marker[sub_type_tag] ,\
                           ls ='--')               

            

            plt.ylim([1.0,100.0])
            plt.xlim([10**10,10**12])
            xlabel('$M_{star} [M_{sun}]$', fontsize= fontsize_plot)

            ylabel('$\Delta \Sigma _{norm}, M^{sub}_{norm} $' , fontsize= fontsize_plot)

            title(sub_type_tag )
            legend(loc='upper left',fontsize=14)


            plot_name = savename+ '_ESDvsMass_' + sub_type_tag 

            plt.savefig(plot_name + ".pdf")
            #plt.savefig(plot_name + ".ps")               
            print 'evince ' + plot_name+ ".pdf"







            

            plt.figure(figsize=(7,6))
            for ibin in np.arange(1 , len(ESD_data.mass_arr)-1) :


                #print i
                # if there are too many error bars reduce the plotted number of pints


                #print ESD_data_list[i].ESD_profiles[sub_type_tag][ibin,plot_index]
                if sub_type_tag == 'sat' :
                    plt.loglog(ESD_data.ESD_profiles['rad'][ibin,:], \
                                       ESD_data.ESD_profiles['cen'][ibin,:] , \
                                       color ='#BDBDBD', \
                                       label='_nolegend_', \
                                       marker=ESD_data.marker[sub_type_tag] ,\
                                       ls =ESD_data.linestyle[sub_type_tag])
                    plt.axvline(x=ESD_data.ESD_profiles['rad'][ibin,index_500kpc], ymin=0., ymax = 1.0, linewidth=2, color=color_cm(0), ls ='--')
                    
                plt.loglog(ESD_data.ESD_profiles['rad'][ibin,:], \
                           ESD_data.ESD_profiles[sub_type_tag][ibin,:] , \
                           color =color_cm(scale_color*ibin), \
                           label=  '$  ' + str(ESD_data.mass_arr[ibin]) + "<\mathrm{log}_{10}\,M_{\mathrm{star}}\, / \, M_{\\odot}<"  + str(ESD_data.mass_arr[ibin+1])+'$  ',\
                           marker=ESD_data.marker[sub_type_tag] ,\
                           ls =ESD_data.linestyle[sub_type_tag])
                
                plt.axvline(x=ESD_data.ESD_profiles['rad'][ibin,index_50kpc], ymin=0., ymax = 1.0, linewidth=2, color=color_cm(0), ls ='--')



            plt.ylim([1.0,1000.0])
            plt.xlim([0.01,2.0])
            xlabel('$R \, [h^{-1}Mpc]$', fontsize= fontsize_plot)

            ylabel('$\Delta \Sigma \, [h M_{\\odot} pc^{-2}]$ ' , fontsize= fontsize_plot)

            if sub_type_tag == 'all':
                title('Centrals + Satellites',fontsize=25)
            elif sub_type_tag == 'cen' :
                title('Centrals',fontsize=25)
            elif sub_type_tag == 'sat' :
                title('Satellites',fontsize=25)  

            legend(loc='upper right',fontsize=12)

            plot_name = savename+ '_allbins_' + sub_type_tag 

            plt.savefig(plot_name + ".pdf")
            #plt.savefig(plot_name + ".ps")               
            print 'evince ' + plot_name+ ".pdf"

            

def compute_chi(ESD_data_list, path):
        sub_type = ['all' , 'sat' ,'cen']
        for ESD_data in ESD_data_list :
                if isinstance(ESD_data, ESD_gama):
                       
                        ESD_gama_comp = ESD_data
        for ESD_data in ESD_data_list :
            for sub_type_tag in sub_type :
                
                if isinstance(ESD_data, ESD_sim):
                    chi_sq_red , dof = ESD_data.compute_chi_squared(ESD_gama_comp,sub_type_tag, path )
                    print   '========== ', ESD_data.name_tag , sub_type_tag ,\
                        'X^2='+ str(chi_sq_red),\
                        'pvalue=' , 1 - stats.chi2.cdf(chi_sq_red*dof,dof ),\
                         chi2.sf(chi_sq_red*dof,dof)
def plot_ESD_list(ESD_data_list, savename):
        """ plot all esd profiles in one plot"""
        # if there are more than one sim ESD use def colors, with one use blue = sat red = cen black = all
        n_ESD_sim = 0
        n_ESD_gama =0
        for ESD_data in ESD_data_list :
                if isinstance(ESD_data, ESD_sim):
                        n_ESD_sim +=1
                        ESD_sim_comp = ESD_data
                if isinstance(ESD_data, ESD_gama):
                        n_ESD_gama +=1
                        ESD_gama_comp = ESD_data

        subplotpick = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]    
                        
        scale_color = 40
        sub_type = ['all' , 'sat' ,'cen']
        
        sub_type_cm = [plt.cm.gray , plt.cm.winter ,plt.cm.autumn]
        #----------- plotting the profiles in diff bins
        rcParams['figure.figsize'] = 20, 15

        #----------------- plotting the ESD profiles
        scale_color = 0
        for sub_type_tag , color_cm  in zip(sub_type , sub_type_cm):
                plt.figure(figsize=(18,11))
                for pbin in np.arange(len(ESD_data_list[0].mass_arr)-2) :
                    ibin = pbin +1
                    panel = plt.subplot2grid((2,3), subplotpick[pbin], \
                    colspan=1,rowspan=1)
            # 		We change the fontsize of minor ticks label 
                    plt.tick_params(axis='both', which='major', labelsize=20)
                    plt.tick_params(axis='both', which='minor', labelsize=20)
                    loglog( nonposy='clip')
                    xlim(0.01,2.1)
                    ylim(2,900)
                    subplots_adjust(hspace=0)
                    subplots_adjust(wspace=0)
                    panel.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
                    plt.xticks((0.10, 1.0), ('0.10', '1.0') )
                    #panel.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
                    if(pbin==0):
                        text(0.003, 2,'$\Delta \Sigma \, \, [h M_{\\odot}/{\\rm pc}^2]$', ha='center', va='center',rotation='90',fontsize=35)
                        panel.set_xticks([])
                    elif(pbin==1)|(pbin==2)|(pbin==5):
                        panel.set_yticks([])
	            elif(pbin==4):
		        text(0.2  ,1,'R $\, [h^{-1} {\\rm Mpc}]$', ha='center', va='center',fontsize=35)
                        panel.set_yticks([])

                

        

                    for i in range(0, size(ESD_data_list)) :
                            #print i
                        # if there are too many error bars reduce the plotted number of pints
                            if n_ESD_sim == 1 and n_ESD_gama == 1: color_line =color_cm(scale_color*ibin)
                            else: color_line =ESD_data_list[i].color[sub_type_tag]
                            n_points = size(ESD_data_list[i].ESD_profiles['rad'][ibin,:])
                            if (n_points  >  20.0) :

                                    plot_index = range(i , n_points , 5)

                            else:
                                    plot_index = range(0 , n_points )

                            #print ESD_data_list[i].ESD_profiles[sub_type_tag][ibin,plot_index]
                            errorbar(ESD_data_list[i].ESD_profiles['rad'][ibin,plot_index], \
                                         ESD_data_list[i].ESD_profiles[sub_type_tag][ibin,plot_index] , \
                                         yerr= [(ESD_data_list[i].ESD_profiles[sub_type_tag][ibin,plot_index] - ESD_data_list[i].ESD_16th[sub_type_tag][ibin,plot_index]) \
                                                , (ESD_data_list[i].ESD_84th[sub_type_tag][ibin,plot_index] - ESD_data_list[i].ESD_profiles[sub_type_tag][ibin,plot_index])] , \
                                         color =color_line, \
                                         label=ESD_data_list[i].name_tag ,\
                                         marker=ESD_data_list[i].marker[sub_type_tag] ,\
                                         ls =ESD_data_list[i].linestyle[sub_type_tag])


                                

                    text(0.2 , 600 , '$  ' + str(ESD_data_list[0].mass_arr[ibin]) + "<\mathrm{log}_{10}\,M_{\mathrm{star}}\, / \, M_{\\odot}<"  + str(ESD_data_list[0].mass_arr[ibin+1])+'$  ', ha='center', va='center',fontsize=20)

                    if pbin == 4 :legend(loc='center left',fontsize=18, bbox_to_anchor=(1, 0.5))

                if sub_type_tag == 'all':
                    plt.suptitle('Centrals + Satellites',fontsize=30)
                elif sub_type_tag == 'cen' :
                    plt.suptitle('Centrals',fontsize=30)
                elif sub_type_tag == 'sat' :
                    plt.suptitle('Satellites',fontsize=30)                
                    
                plot_name = savename+  '_multi_' + sub_type_tag 

                plt.savefig(plot_name + ".pdf")
                #plt.savefig(plot_name + ".ps")

                print 'evince ' + plot_name+ ".pdf"







        #----------------- plotting the F_sat

        plt.figure(figsize=(7,6))

        for i in range(size(ESD_data_list)) :
                    plt.plot(ESD_data_list[i].mean_mass[1:] , ESD_data_list[i].f_sat[1:], color =ESD_data_list[i].color['all'], label=ESD_data_list[i].name_tag , marker=ESD_data_list[i].marker['all'] , ls =ESD_data_list[i].linestyle[sub_type_tag])

        ylim(0,1)

        xlabel('$\mathrm{log}_{10} \, M_{star} \, / \, M_{\\odot}$', fontsize= fontsize_plot)

        ylabel('$f_{sat}$ ' , fontsize= fontsize_plot)

        #title(sub_type_tag + '  ' +min_mass_tag + "<M<" + max_mass_tag )

        legend(loc='lower left',fontsize=14)#'upper right')

        plot_name = savename+'_f_sat' 
        plt.savefig(plot_name+ '.pdf')
        #plt.savefig(plot_name+ '.ps')        
        print 'evince ' + plot_name

        


def plot_ESD_sim_stat(ESD_data_list, savename, save_txt = False):
        """ plots the statistics in every bin"""
        # if there are more than one sim ESD use def colors, with one use blue = sat red = cen black = all
        n_ESD_sim = 0
        n_ESD_gama =0
        ESD_sim_comp = []
        ESD_gama_comp = []
        for ESD_data in ESD_data_list :
                if isinstance(ESD_data, ESD_sim):
                        n_ESD_sim +=1
                        ESD_sim_comp.append(ESD_data)
                if isinstance(ESD_data, ESD_gama):
                        n_ESD_gama +=1
                        ESD_gama_comp.append(ESD_data)


        

        
                        
        for mass_type_tag in ESD_sim_comp[0].avail_stat:

        #----------------- plotting the ESD profiles

            for sub_type_tag  in ESD_sim_comp[0].sub_type :
                    plt.figure(figsize=(18,11))
                    for ibin in np.arange(len(ESD_sim_comp[0].mass_arr)-1) :
                            ax1 = plt.subplot(2, 3, ibin+1)#, sharex=True, sharey=True)

                            min_x = np.array([])
                            max_x =np.array([])

                            for ESD_s in ESD_sim_comp :

                                
                                if save_txt :
                                    filename =savename+'_'+ ESD_s.name +'_'+ sub_type_tag+'_'+ mass_type_tag+ str(ESD_s.mass_arr[ibin]) + "_"  + str(ESD_s.mass_arr[ibin+1])+'.txt'
                                    print 'less ' , filename 
                                    with open(filename, "w") as text_file:
                                        text_file.write(("middlebin , "+\
                                                         "density \n"))
                                        print len(ESD_s.vec_stat[ibin][sub_type_tag]['histo_x_'+ mass_type_tag][0])
                                        for bin_histo in np.arange(len(ESD_s.vec_stat[ibin][sub_type_tag]['histo_x_'+ mass_type_tag][0])-1) :
                                            
                                            text_file.write((" {:03.2f} ,"+\
                                                             " {:03.2f} \n")\
                                                            .format(ESD_s.vec_stat[ibin][sub_type_tag]['histo_x_'+ mass_type_tag][0,bin_histo],\
                                                                    ESD_s.vec_stat[ibin][sub_type_tag]['histo_y_'+ mass_type_tag][0,bin_histo]))

                                
                                x_arr = ESD_s.vec_stat[ibin][sub_type_tag]['histo_x_'+mass_type_tag][0,:]
                                if min(x_arr) != max(x_arr):
                                    min_x = np.append(min_x , min(x_arr))
                                    max_x = np.append(max_x , max(x_arr))            
                                
                                    ax1.plot(ESD_s.vec_stat[ibin][sub_type_tag]['histo_x_'+mass_type_tag][0,:], \
                                             ESD_s.vec_stat[ibin][sub_type_tag]['histo_y_'+mass_type_tag][0,:] , \
                                             color =ESD_s.color[sub_type_tag],\
                                             label=ESD_s.name_tag + 'N = ' + str(ESD_s.vec_stat[ibin][sub_type_tag]['n_sub']),\
                                             marker=ESD_s.marker[sub_type_tag] ,\
                                             ls =ESD_s.linestyle[sub_type_tag])
                                    ax1.axvline(x=ESD_s.vec_stat[ibin][sub_type_tag]['median_'+mass_type_tag], ymin=0., ymax = 1.0, linewidth=2, color=ESD_s.color[sub_type_tag], ls ='--')

                            #ax1.set_yscale('log', nonposy='clip')
                            #ax1.set_xscale('log')

                            
                            #plt.ylim([1.0,1000.0])
                            plt.xlim([ min(min_x) , max(max_x) ])
                            
                            #xlabel('$R \, [h^{-1}Mpc]$', fontsize= fontsize_plot)

                            #ylabel('$\Delta \Sigma \, [h M_{sun} pc^{-2}]$ ' , fontsize= fontsize_plot)

                            title(sub_type_tag + '  ' + str(ESD_sim_comp[0].mass_arr[ibin]) + "<M<"  + str(ESD_sim_comp[0].mass_arr[ibin+1]))

                            legend(loc='upper right')

                    plt.suptitle(mass_type_tag)
                    plot_name = savename+  '_stat_'+ mass_type_tag+'_of_' + sub_type_tag  

                    plt.savefig(plot_name + ".pdf")
                    #plt.savefig(plot_name + ".ps")
                    #plt.show()
                    print 'evince ' + plot_name+ ".pdf"







                    


def print_stat_table(ESD_data_list, savename):
        """print the table with the important stat in latex format"""
        #check that is a sim data pair
        #print only values
        n_ESD_sim = 0
        n_ESD_gama =0
        for ESD_data in ESD_data_list :
                if isinstance(ESD_data, ESD_sim):
                        n_ESD_sim +=1
                        ESD_sim_vec = ESD_data
                if isinstance(ESD_data, ESD_gama):
                        n_ESD_gama +=1
                        ESD_gama_vec= ESD_data
        if (n_ESD_sim != 1) | (n_ESD_gama != 1) : 
            raise NameError('there is more than a couple of data end obs.. exiting')
        with open(savename+".txt", "w") as text_file:
            text_file.write(("$M_{\rm star}$ & "+\
                             " $  M_{200}^{\rm crit}|_{\rm cen}$ & "+\
                             "$\  M_{200}^{\rm crit}|_{\rm sat}$ & "+\
                             "$\ M_{\rm sub}^{\rm cen}$ & $  "+\
                             "M_{\rm sub}^{\rm sat}$ & "+\
                             "$ M_{\rm sub}^{\rm sat}/M_{200}^{\rm crit}$ & "+\
                             "$ d_{\rm sat}$ & $  r_{\rm half}^{\rm dm}|_{\rm cen}$ & "+\
                             "$  r_{\rm half}^{\rm dm}|_{\rm sat}$ & "+\
                             "$N_{\rm gal}$ & "+\
                             "$  M_{\rm lim}^{\rm EAGLE}$ & "+\
                             "$  f_{\rm sat}^{}$   \\\\ \n"))
            
            text_file.write(("* & ** & ** & ** & ** &  & *** & *** & *** &   & **& \\\\ \n"))
            text_file.write(("(1)   &(2)            & (3)     & (4)    & (5)    & (6)    & (7)    &(8)  &(9)  &(10)&(11)&(12)\\\\ \n"))
            for ibin in np.arange(len(ESD_sim_vec.mass_arr)-1) :
                text_file.write(("$[{}-{}] $ & "+\
                                 "$ {:03.2f} $ & "+\
                                 "$ {:03.2f} $ & "+\
                                 "$ {:03.2f} $ & "+\
                                 "$ {:03.2f} $ & "+\
                                 "$ {:03.2f} $ & "+\
                                 "$ {:03.1f} $ & "+\
                                 "$ {:03.1f} $ & "+\
                                 "$ {:03.1f} $ & "+\
                                 "$ {} $ & "+\
                                 "$ {:03.2f} $ & "+\
                                 "$ {:03.2f} $  "+\
                                 "\\\\  \n")\
                                .format(ESD_sim_vec.mass_arr[ibin] , ESD_sim_vec.mass_arr[ibin+1] , \
                                        ESD_sim_vec.vec_stat[ibin]['cen']['mean_main_halo_m200_crit'][0],\
                                        ESD_sim_vec.vec_stat[ibin]['sat']['mean_main_halo_m200_crit'][0],\
                                        ESD_sim_vec.vec_stat[ibin]['cen']['mean_subhalo_mass'][0],\
                                        ESD_sim_vec.vec_stat[ibin]['sat']['mean_subhalo_mass'][0],\
                                        10**ESD_sim_vec.vec_stat[ibin]['sat']['mean_Msub_Mhost200_ratio'][0],\
                                        
                                        10**(ESD_sim_vec.vec_stat[ibin]['sat']['mean_sub_halo_rad_dist_main'][0]+3),\
                                        10**(ESD_sim_vec.vec_stat[ibin]['cen']['mean_subhalo_half_mass_rad_dm'][0]+3),\
                                        10**(ESD_sim_vec.vec_stat[ibin]['sat']['mean_subhalo_half_mass_rad_dm'][0]+3),\
                                        ESD_sim_vec.vec_stat[ibin]['all']['n_sub'][0],\
                                        ESD_sim_vec.m_host_limit[ibin],\
                                        float(ESD_sim_vec.vec_stat[ibin]['sat']['n_sub'][0])/ float(ESD_sim_vec.vec_stat[ibin]['all']['n_sub'][0]),\
                                        #ESD_gama_vec.f_sat[ibin]
                                        
                ))

        print 'less ' , savename+".txt"
        
                        




def save_stat_table(ESD_data_list, savename):

        """ save a table with the median and percentile stats """
        #print the table with the important stat in latex format
        #check that is a sim data pair
        #print only values
        n_ESD_sim = 0
        n_ESD_gama =0
        for ESD_data in ESD_data_list :
                if isinstance(ESD_data, ESD_sim):
                        n_ESD_sim +=1
                        ESD_sim_vec = ESD_data
                if isinstance(ESD_data, ESD_gama):
                        n_ESD_gama +=1
                        ESD_gama_vec= ESD_data
        if (n_ESD_sim != 1) | (n_ESD_gama != 1) : 
            raise NameError('there is more than a couple of data end obs.. exiting')
        with open(savename+"_marcello.txt", "w") as text_file:
            text_file.write(("Mstar_median , "+\
                             "Mstar_16p , "+\
                             "Mstar_84p , "+\
                             "Msub_median , "+\
                             "Msub_16p , "+\
                             "Msub_84p  \n"))
            
            for ibin in np.arange(len(ESD_sim_vec.mass_arr)-1) :
                text_file.write((" {:03.2f} ,"+\
                                 " {:03.2f} ,"+\
                                 " {:03.2f} ,"+\
                                 " {:03.2f} ,"+\
                                 " {:03.2f} ,"+\
                                 " {:03.1f} \n")\
                                .format(ESD_sim_vec.vec_stat[ibin]['all']['median_subhalo_mass_in_star'][0],\
                                        ESD_sim_vec.vec_stat[ibin]['all']['16th_p_subhalo_mass_in_star'][0],\
                                        ESD_sim_vec.vec_stat[ibin]['all']['84th_p_subhalo_mass_in_star'][0],\
                                        ESD_sim_vec.vec_stat[ibin]['all']['median_subhalo_mass'][0],\
                                        ESD_sim_vec.vec_stat[ibin]['all']['16th_p_subhalo_mass'][0],\
                                        ESD_sim_vec.vec_stat[ibin]['all']['84th_p_subhalo_mass'][0],\
                                        ))
  

        print 'less ' , savename+"_marcello.txt"


        
    

import copy as cp
import pyfits
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import chisquare
from scipy import stats
from scipy.stats.distributions import chi2

#color_list
gray_1 = "#a8a8a8"
gray_2 = "#7e7e7e"
gray_3 = "#939393"




#mpl.rcParams['font.size']       = 14.0
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['axes.linewidth']  =  1.5

fontsize_plot =20



#====================== Paths ++++++++++++++++++

#path where to put the plots
path_plot = "/disks/galform11/velliscig/Documents/Draft4/py_figs/"
#create a list that contains all the objects to plot


#path for the catalogue from the simulations
path_sim_z0p2 = "/net/galaxy/data1/velliscig/project2/REFERENCE_PLANCK_L0100N1504/z0p18/data/REFERENCE_PLANCK_L0100N1504_ESD_stat_Z_MsunM_host_limit_"


#path where the KIDS ESD data is
path_gama_data = '/net/zoom/data2/dvornik/test/test_marcello/COV_ESD_EAGLE_ALL/output_logmstar_bins_Z_B0p005-1p2_Rbins10-20-2000kpc_Om0p315_Ol0p685_Ok0_h1/results_shearcovariance/'

# extra tag needed for the KIDS filename
extra_tag = 'IDs12248'


# vector that specify the binning of the mass range
# this is in log10(M_star [Msun])
#!!!! this vectohas to correspond to the one used for the KIDS pipeline calculation
massbin_vector =  [10.,10.3,10.6,10.9,11.2,11.5,11.8]




#====================== Calculations  ++++++++++++++++++

# vector of objects that passed to the plotting routine
ESD_data_list_allstars = []

#multiplicative bias factor for the ESD profile correction
m_bias=-0.012

orange = "#E78532"
# initialization of the object ESD_gama
#  ESD_gama( mass_arr , path, nfof_lim , color, name_tag, extra_tag , cat_blind , m_bias)
data_g= ESD_gama('' ,massbin_vector , path_gama_data  , 5 , 'purple', "Blind 2", extra_tag , 'C', m_bias)
data_g.redefine_markers({'all': '^'  , 'cen': 'd'  , 'sat': 's' })
ESD_data_list_allstars.append(data_g)



data_g= ESD_gama('' ,massbin_vector , path_gama_data  , 5 , "green", "Blind 3", extra_tag , 'A' , m_bias)
data_g.redefine_markers({'all': '^'  , 'cen': 'd'  , 'sat': 's' })
ESD_data_list_allstars.append(data_g)

data_g= ESD_gama('' ,massbin_vector , path_gama_data  , 5 , "k", "Blind 1", extra_tag , 'B', m_bias)
data_g.redefine_markers({'all': '^'  , 'cen': 'd'  , 'sat': 's' })
ESD_data_list_allstars.append(data_g)


calib_vect_Mlim = np.array( [9.33867    ,  9.46250  ,    9.91261 ,     9.95804   ,   10.3300   ,   9.5000])


# ESD_sim(mass_arr, path , nfof_lim, color, name_tag, m_host_limit, extra_tag)
data_s = ESD_sim('calibrated_Mlim' ,massbin_vector ,path_sim_z0p2 , 5 , orange , "EAGLE" ,  calib_vect_Mlim , '')
data_s.redefine_colors({'all': 'k'  , 'cen': 'r'  , 'sat': 'b' })
data_s.read_EAGLE_gal_cat()
data_s.stat_on_bins()
data_s.ESD_in_bins()

ESD_data_list_allstars.append(data_s)





for ESD_class in ESD_data_list_allstars:

        if isinstance(ESD_class, ESD_sim):
            # for the simulations read the catalogue
            ESD_class.read_EAGLE_gal_cat()
            # compute profiles and statistics in bins
            ESD_class.stat_on_bins()
        # read f_sat from files or compute it (sim)
        ESD_class.get_fsat()
        # read the ESD profile from file or compute it (sim)
        ESD_class.get_ESD_profile()


savename = path_plot + "allstars"

plot_ESD_sim_stat(ESD_data_list_allstars, savename,save_txt=True)

# plots with all massbins profiles in one plot (only for sim)
plot_ESD_allbins(data_s, savename)
# plots with the comparison data sim one panel for every mass bin
plot_ESD_list(ESD_data_list_allstars, savename)





#make the table in latex for paper
#print_stat_table(ESD_data_list_allstars, savename )
# make text files with the stats
#save_stat_table(ESD_data_list_allstars, savename )

#check problem with errors
ESD_data_list_allstars[3].ESD_prof_in_bin[0]['sat']['boot_84th_snap_esd_profile_y'] - ESD_data_list_allstars[3].ESD_prof_in_bin[0]['sat']['boot_16th_snap_esd_profile_y'] 
ESD_data_list_allstars[3].ESD_84th[0]['sat'] -  ESD_data_list_allstars[3].ESD_16th[0]['sat']
#!!!!!!!!!!!!!!!!!!!!! example of use stat in bin !!!!!!!!!!!!!
#histogram of the distribution of subhalo masses for satellites in the first bin
#data_s.vec_stat[0]['sat']['histo_y_subhalo_mass']
#data_s.vec_stat[0]['sat']['histo_x_subhalo_mass']






# #---------------------------------------------------------------------------------------------------------
print 'new'
#plot_ESD_list(ESD_data_list_diff_Mlim, savename)
#stop()
ESD_data_list_diff_Mlim = []





data_s = ESD_sim('calibrated_Mlim' , massbin_vector ,path_sim_z0p2 , 5 , orange , "EAGLE" ,  calib_vect_Mlim , '')


data_s = ESD_sim('calibrated_Mlim_minus0p75' , massbin_vector ,path_sim_z0p2 , 5 ,    'b', '$\mathrm{log}_{10}\,M^{\mathrm{limit}}_{\mathrm{fiducial}} -0.75 $' , calib_vect_Mlim[:] -0.75, '' )
data_s.redefine_colors({'all': gray_3  , 'cen': gray_3  , 'sat': gray_3 })
#data_s.redefine_linestyles({'all': '--'  , 'cen': '--'  , 'sat': '--' })
#data_s.redefine_linestyles(self,dic)


ESD_data_list_diff_Mlim.append(data_s)

data_s = ESD_sim('calibrated_Mlim_plus0p25' , massbin_vector ,path_sim_z0p2 , 5 ,    'r', '$\mathrm{log}_{10}\,M^{\mathrm{limit}}_{\mathrm{fiducial}} + 0.25 $' , calib_vect_Mlim +0.25, '' )
data_s.redefine_colors({'all': gray_1  , 'cen': gray_1  , 'sat': gray_1 })
#data_s.redefine_linestyles({'all': '-.'  , 'cen': '-.'  , 'sat': '-.' })

ESD_data_list_diff_Mlim.append(data_s)

#data_s = ESD_sim(path_sim_z0p2 , 4 ,    blue_alt, " Calib -1.0 " , calib_vect_Mlim -1.0, '' )
#ESD_data_list_diff_Mlim.append(data_s)

data_s = ESD_sim('calibrated_Mlim_minus1p50' , massbin_vector ,path_sim_z0p2 , 5 ,    'g',  '$\mathrm{log}_{10}\,M^{\mathrm{limit}}_{\mathrm{fiducial}} -1.5 $' , calib_vect_Mlim -1.5, '' )
data_s.redefine_colors({'all': gray_2  , 'cen':  gray_2 , 'sat': gray_2 })
#data_s.redefine_linestyles({'all': ':'  , 'cen': ':'  , 'sat': ':' })
ESD_data_list_diff_Mlim.append(data_s)

data_s = ESD_sim('calibrated_Mlim' , massbin_vector ,path_sim_z0p2 , 5 ,    'k', '$\mathrm{log}_{10}\,M^{\mathrm{limit}}_{\mathrm{fiducial}}$', calib_vect_Mlim  , '' )
data_s.redefine_colors({'all': 'k'  , 'cen': 'r'  , 'sat': 'b' })
ESD_data_list_diff_Mlim.append(data_s)


data_g= ESD_gama('' , massbin_vector , path_gama_data  , 5 , "k", "KiDS $N^{\mathrm{GAMA}}_{\mathrm{FoF}}	\geqslant 5$", extra_tag , 'B', m_bias)

#data_g= ESD_gama(path_gama + "cat_v2_vol_limit_fluxcscale_corr/" , 4, 'k', "KiDS nFOF>4")
#data_g.redefine_colors({'all': 'k'  , 'cen': 'r'  , 'sat': 'b' })
data_g.redefine_markers({'all': '^'  , 'cen': 'd'  , 'sat': 's' })
ESD_data_list_diff_Mlim.append(data_g)



for ESD_class in ESD_data_list_diff_Mlim:

        if isinstance(ESD_class, ESD_sim):
            ESD_class.read_EAGLE_gal_cat()
            ESD_class.stat_on_bins()
            ESD_class.ESD_in_bins()
            
        ESD_class.get_fsat()
        ESD_class.get_ESD_profile()
        



                                               
savename = path_plot + "diff_Mlim"
plot_ESD_list(ESD_data_list_diff_Mlim, savename)
#plot_ESD_sim_stat(ESD_data_list_diff_Mlim, savename)
compute_chi(ESD_data_list_diff_Mlim, path_plot)

plot_ESD_list(ESD_data_list_diff_Mlim, savename)


stop()















# #---------------------------------------------------------------------------------------------------------


ESD_data_list_allstars = []


#define the objects ESD_gama and append them in a list


data_g= ESD_gama(path_gama + "cat_v2_vol_limit_fluxcscale_corr/" , 0, "red", "Mstar fluxscale nFOF >1")
data_g.get_fsat()

ESD_data_list_allstars.append(data_g)


orange = "#E78532"
data_s = ESD_sim('' ,path_sim_z0p2 , 0 , orange , "all stars nFOF >1" ,  [ 10.2239  ,    10.5137   ,   10.8012 ,     11.0813    ,  11.3906 ,     8.00000] , '')
data_s.read_EAGLE_gal_cat()
data_s.stat_on_bins()
data_s.ESD_in_bins()
#data_s = ESD_sim(path_sim_z0p2 , 4 , orange , "all stars nFOF >1" ,   [ 10.,  10.   ,  10.   ,   10.    ,  10.   ,   10.] , '')
# data_s = ESD_sim(path_sim_z0p2 , 4 , orange , "all stars nFOF >1" ,   [ 8.5,  8.5   ,  8.5   ,   8.5    ,  8.5   ,   8.5] , '')
# data_s = ESD_sim(path_sim_z0p2 , 4 , orange , "all stars nFOF >1" ,   [ 8.,  8.   ,  8.   ,   8.    ,  8.   ,   8.] , '')

data_s.read_EAGLE_gal_cat()

data_s.ESD_in_bins()



ESD_data_list_allstars.append(data_s)

for ESD_class in ESD_data_list_allstars:
                                      
        ESD_class.get_fsat()
        ESD_class.get_ESD_profile()
                                      

#print data_s.compute_chi_squared(data_g, 'all')

#print data_s.ESD_prof_in_bin[3]['all']['mean_snap_esd_profile_y']
#print data_s.ESD_prof_in_bin[3]['all']['snap_r_plus_half']

#print data_s.ESD_profiles[3]['all']
#print data_s.ESD_profiles[3]['rad']

#print abs(data_s.ESD_profiles[3]['all'] - data_s.ESD_prof_in_bin[3]['all']['mean_snap_esd_profile_y'])
#print abs(data_s.ESD_16th_error[3]['all'] + data_s.ESD_profiles[3]['all'] - data_s.ESD_prof_in_bin[3]['all']['boot_05th_snap_esd_profile_y']/100)


#quit()       
savename = path_plot + "allstars_nFOF0"
plot_ESD_list(ESD_data_list_allstars, savename)


# #---------------------------------------------------------------------------------------------------------

ESD_data_list_30Kpc = []
        
data_g3 = ESD_gama(path_gama + "cat_v2_vol_limit_fluxscale_30kpc" , 4, 'blue' , "Mstar 30kpc corr")

ESD_data_list_30Kpc.append(data_g3)
                                                                      
blue_alt = "#4477AA"
data_s4 = ESD_sim(path_sim_z0p2  , 4,   blue_alt, "30Kpc calib " , [        9.34   ,   9.58    ,  10.03    ,  10.37   ,   10.35  ,    10] , 'Mstar30kpc' )
ESD_data_list_30Kpc.append(data_s4)

for ESD_class in ESD_data_list_30Kpc:
        ESD_class.get_fsat()
        ESD_class.get_ESD_profile()

savename = path_plot + "30kpc"
plot_ESD_list(ESD_data_list_30Kpc, savename)
compute_chi(ESD_data_list_30Kpc)

# #---------------------------------------------------------------------------------------------------------


ESD_data_list_diff_Nfof = []

data_g= ESD_gama(path_gama + "cat_v2_vol_limit_fluxcscale_corr/" , 4, "red", "Mstar fluxscale nFOF >4")

ESD_data_list_diff_Nfof.append(data_g)

data_s = ESD_sim(path_sim_z0p2 , 0 ,    'k', " Nfof>0 " , [9.0 , 9.0 , 9.0 , 9.0 , 9.0 , 9.0] , '' )
ESD_data_list_diff_Nfof.append(data_s)
data_s = ESD_sim(path_sim_z0p2 , 2 ,    'b', " Nfof>2 " , [9.0 , 9.0 , 9.0 , 9.0 , 9.0 , 9.0] , '' )
ESD_data_list_diff_Nfof.append(data_s)
data_s = ESD_sim(path_sim_z0p2 , 4 ,    'r', " Nfof>4 " , [9.0 , 9.0 , 9.0 , 9.0 , 9.0 , 9.0] , '' )
ESD_data_list_diff_Nfof.append(data_s)
data_s = ESD_sim(path_sim_z0p2 , 6 ,    'g', " Nfof>6 " , [9.0 , 9.0 , 9.0 , 9.0 , 9.0 , 9.0] , '' )
ESD_data_list_diff_Nfof.append(data_s)
data_s = ESD_sim(path_sim_z0p2 , 8 ,    'y', " Nfof>8 " , [9.0 , 9.0 , 9.0 , 9.0 , 9.0 , 9.0] , '' )
ESD_data_list_diff_Nfof.append(data_s)

for ESD_class in ESD_data_list_diff_Nfof:
        ESD_class.get_fsat()
        ESD_class.get_ESD_profile()
        if isinstance(ESD_class, ESD_sim):
            ESD_class.read_EAGLE_gal_cat()
            ESD_class.stat_on_bins()        


        
savename = path_plot + "diff_Nfof"
plot_ESD_list(ESD_data_list_diff_Nfof, savename)
plot_ESD_sim_stat(ESD_data_list_diff_Nfof, savename)

















