#MAKEFLAGS += -j5

### ----------------------------------
# This file create directional waves spectra from ICEsat2 tracks.
# The algorithm is devided in 3 subsections and this nomencalture
# is also used to defined input and output file names, and figures.
#
# How to use this file:
# This file can operate on batches of track IDs. Once a batch is created
# from a list of track IDs (using the "downlaod"), all susequent targets
# can be generated. The tracklist can be modified by just linked to a new track_list.
# the processing steps are in order, such that higher numbers generally depend on (most) of the lower numbers.
#
# Dependencies:
# This uses the XX libary () and .netrc authentification ot communicate wit the NSIDC server.
#
# A) download data:
#
# A01a Downlaod ALT07 tracks given a track ID list or a single track ID:
#    usage with track ID list:
#    make -f make.make batch=batch_name track_list=track_ID_list.json A01a
#    track_ID_list.json must be in folder ../track_lists/
#    use for single ALT03 track_ID
#    make -f make.make batch=batch_name track_ID=track_ID A01a
#    or, if existing file should be used;
#    make -f make.make batch=batch_name A01a
#    this corresponds to track_ID=False and track_list=False
#
# A01b test if there can be waves in track, decides to downlaod ALT03 tracks
#    and generates approxiamte filenames
# A01c downlaods corrsponding ATL03 tracks
#
# A02 generate prior data from WW3 data, for a given, existing ATL03 track
# A03 download corrsponding ALT10 track (not done jet)
#
# B) derive spectra and directional estimates
# B01 				filters and regridds data, data is saves in binnned .h5 tables.
# B02 				makes gFT and FFT spectra
# B02_collect collects the plots in a new folder
# B03 				plots non-correctd spectral estimates
# reset_B03		resets B03 targets
# B04					MCMC angle inversion based on A02b
# B05					Filter for best guess incident angle and plot overview
# B06					correct spectral estimates, calculates continues residual, and reformats the data.
# 						all results are storred and collected in one place
#
# C) plotting routines (not implemented yet)
# P) plotting routines to make piublishable figures (not implimented yet)
#
# ### helpers:
# print-var_name  prints the variable in output
# example         shows example configurations for executing the make file
#	make_folder			makes (nearly) all folders needed for the procesing based on the defined paths
# rm_bad_beams    removes beams that are not sufficient from B02 file list (ALT03 tracks)
# sync_gdrive 		syncs results from batch_key to gdrive, needs rclone to be configured
#
#
# options and flags:
# -n dry runs
# -B forces to rebuild
# -jN number of parallel processes executed
# -k, --keep-going
# 		 Continue  as  much  as  possible after an error.  While the target
# 		 that failed, and those that depend on it, cannot  be  remade,  the
# 		 other dependencies of these targets can be processed all the same.
#
# -l [load], --load-average[=load]
# 		 Specifies that no new jobs (commands) should be started  if  there
# 		 are  others  jobs running and the load average is at least load (a
# 		 floating-point number).  With no argument, removes a previous load
# 		 limit.
### ----------------------------------

### example : show example how to use this make file
.PHONY : example
example :
					@echo "examples:"
					@echo " "
					@echo "simple target excecution: "
					@echo "make batch=batch02 track_list=batch02_tracks_test.json A01a -n"
					@echo " "
					@echo "to download data given a batch list; -n is the text flag"
					@echo "make batch=batch02 track_list=batch02_tracks_test.json A01 -n"
					@echo " "
					@echo "to make targets without re-downloading the data; -n is the text flag"
					@echo "make batch=batch02 track_list=False B01 -n"
					@echo " "
					@echo "print targets"
					@echo "make batch=batch02 track_list=False print-B01_targets -n"
					@echo " "
					@echo "force target to rebuld -B"
					@echo "make batch=batch02 track_list=False B01 -B -n"

# define subroutines
MKDIR_P = mkdir -p

# define paths based on config.json file in ../config/
analysisfolder		= $(shell jq -r '.paths.base' ../config/config.json)analysis_db/#/home/mhell/2021_ICESat2_tracks/analysis_db/
plot_folder	     	= $(shell jq -r '.paths.plot' ../config/config.json)#/home/mhell/2021_ICESat2_tracks/plots/
work_folder				= $(shell jq -r '.paths.work' ../config/config.json)#/work/mhell_work/2021_ICESat2_tracks/
scratch_folder    = $(shell jq -r '.paths.scratch' ../config/config.json)#/scratch/mhell/2021_ICESat2_tracks/
track_lists_folder= $(shell jq -r '.paths.base' ../config/config.json)track_lists/#/home/mhell/2021_ICESat2_tracks/track_lists/
track_downloader  = $(shell jq -r '.paths.local_script' ../config/config.json)read-ICESat-2/scripts/#/home/mhell/2021_ICESat2_tracks/modules/read-ICESat-2/scripts/
bad_tracks_folder = $(work_folder)bad_tracks


### ------------ set default parameters -----------------
# This values are used if not defiend otherwise in make command:
# hemisphere or geographical subset of the analysis, currently SH or NH
hemis = SH
# batch
batch = batch_dummy#batch01

# (False) if not False, this should be a track_ID.
# for processing single track_ID ..
track_ID 	= False

# if not False provide a filename.json read this file and download the data. check that you are on the right batch!!
track_list=False

# flag for testing. not used yet.
ID_flag	= True
### ------------- end default parameter --------------------

# infer parameters
batch_key = $(hemis)_$(batch)
# infer hemis number, needed for ALT07/10 tracks
ifeq ($(hemis), SH)
hemis_num = 02
else
hemis_num = 01
endif

# make folders for batch
.PHONY : make_folder
make_folder :
					${MKDIR_P} $(scratch_folder)/$(batch_key)
					${MKDIR_P} $(bad_tracks_folder)/$(batch_key)
					${MKDIR_P} $(bad_tracks_folder)/$(batch_key)_deleted
					${MKDIR_P} $(A01b_path)
					${MKDIR_P} $(A01c_path)

### -----------  section A -----------------
### download data

ifeq ($(track_list), False) # just used the filed in the folder, do not redownload them
A01a_targets= $(shell ls $(scratch_folder)/$(batch_key)/ATL07-$(hemis_num)_*.h5)
A01a_tracks=$(foreach i, $(basename $(A01a_targets)), $(subst ATL07-$(hemis_num)_,,$(i))  )

# possibilty to tough downloaded files such that everything after is rebuild
# only executed when forced build is used
$(A01a_targets) : $(scratch_folder)/$(batch_key)/ATL07-$(hemis_num)_%.h5 :
					touch $(scratch_folder)/$(batch_key)/ATL07-$(hemis_num)_$*.h5

else
# create target list from from file in track list folder. The target downloads the respective keys from the interweb.
A01a_track_names=$(shell jq -r .[] $(track_lists_folder)$(track_list))
A01a_tracks =  $(foreach i, $(A01a_track_names) , $(subst ATL07-$(hemis_num)_,,$(i))  )
A01a_targets=  $(foreach i, $(A01a_track_names) , $(addprefix $(scratch_folder)/$(batch_key)/, $(i).h5 ) )

# downloads each file from from NISDC
$(A01a_targets) : $(scratch_folder)/$(batch_key)/ATL07-$(hemis_num)_%.h5 :
					python $(track_downloader)/nsidc_icesat2_associated2.py --user mhell@ucsd.edu --netrc ~/.netrc --product ATL07 --directory $(scratch_folder)/$(batch_key)/ -F ATL07-$(hemis_num)_$*.h5

endif

A01b_path=$(work_folder)/$(batch_key)/A01b_ID
A01b_targets= $(foreach i, $(A01a_tracks) ,$(addprefix $(A01b_path)/, A01b_$(i)_success.json ) )

# checks weather or not if its is worth it to downlaod data and then downloads the corresponding ATL03 data
$(A01b_targets) : $(A01b_path)/A01b_%_success.json : $(scratch_folder)/$(batch_key)/ATL07-$(hemis_num)_%.h5 #$(analysisfolder)/A01b_ALT07_SHNH_variance_tester.py
					python $(analysisfolder)/A01b_ALT07_SHNH_variance_tester.py $* $(batch_key) $(ID_flag) > log/A02/$*.txt 2>&1


A01b_success= $(shell ls $(A01b_path)/A01b_ID_*.json)
A01b_ATL03_files= $(foreach i, $(A01b_success), $(shell jq -r '.tracks.ATL03[]' $(i) ) )
A01b_track_list= $(foreach i, $(basename $(A01b_success)), $(subst $(A01b_path)/A01b_ID_,,$(i))  )
#A01b_track_list=  $(foreach i, $(A01b_track_list1), $(subst _005_02,_005_01,$(i))  )

# if beam is not provieded track list must be provided
ifeq ($(track_ID), False)
ID_track_list= $(A01b_track_list)
else
ID_track_list= $(track_ID)
endif

# merge files and save
#A01c_download := $(scratch_folder)/$(batch_key)
A01c_path := $(scratch_folder)/$(batch_key)#$(work_folder)/$(batch_key)/A01c_getATL03_$(hemis)
A01c_targets := $(foreach i, $(ID_track_list) ,$(addprefix $(A01c_path)/, A01c_ATL03_$(i)_corrected.h5 ) )

$(A01c_targets) : $(A01c_path)/A01c_ATL03_%_corrected.h5 : $(A01b_path)/A01b_ID_%.json
					python $(analysisfolder)/A01c_mergefiles.py $* $(batch_key) $(ID_flag) > log/A02/$*.txt 2>&1

#$(scratch_folder)/$(batch_key)/ATL03_%.h5
#cannot link this to $(A01b_path)/A01b_success_%.json because the tracknames differ (the region definiiton changes)
#python $(track_downloader)/nsidc_icesat2_associated.py --user mhell@ucsd.edu --netrc ~/.netrc --product ATL03 --directory $(A01c_download) -F processed_ATL10_$*.h5 -R 2


# @echo $(foreach i, $(t_list), ATL03_$(word 3, $(subst _, ,${i}))_$(word 4, $(subst _, ,${i))*.h5 )

#@echo $(foreach i, $(t_list), ATL03_$(word 3, $(subst _, ,${i}))_$(word 4, $(subst _, ,${i))*.h5 )
#@echo $(word 2, $(subst _, ,$*))*$(word 3, $(subst _, ,$*))*.h5
#@echo $(shell ls $(A01c_download)/$(word 2, $(subst _, ,$*))*$(word 3, $(subst _, ,$*))*.h5)

# $(A01c_targets) : $(A01c_path)/A01c_success_%.txt : #$(A01b_path)/A01b_%.json
# 					$(eval t_list := $(foreach i, $(shell jq -r '.tracks.ATL03_guess[]' $(A01b_path)/A01b_ID_$*.json), processed_$(i).h5 ))
# 					$(eval dummy_list := $(foreach i, $(shell jq -r '.tracks.ATL03_dummy[]' $(A01b_path)/A01b_ID_$*.json),$(i) ))
# 					@echo $(t_list)
# 					@echo $(dummy_list)
#
# 					python $(track_downloader)/nsidc_icesat2_associated.py --user mhell@ucsd.edu --netrc ~/.netrc --product ATL03 --directory $(A01c_download) -F $(t_list) -R 2 -P 2
# 					@echo $(foreach i, $(dummy_list), $(A01c_download)/ATL03_$(word 1, $(subst _, ,$(i)  ))*_$(word 2, $(subst _, ,$(i) ))*.h5 )
# 					@echo $(shell ls $(foreach i, $(dummy_list), $(A01c_download)/ATL03_$(word 1, $(subst _, ,$(i)  ))*_$(word 2, $(subst _, ,$(i) ))*.h5 ) )

#@echo "$(t_list)" > $(A01c_path)/A01c_success_$*.txt

.PHONY : A01a A01b A01c
A01a : make_A01_folders $(A01a_targets)
A01b : $(A01b_targets)
A01c : $(A01c_targets)


make_A01_folders:
					${MKDIR_P} $(A01b_path)

## generate prior data from WW3 data, for a given, existing ATL03 track
A02_targets := $(foreach i, $(ID_track_list) ,$(addprefix $(work_folder)/$(batch_key)/A02_prior/, A02_${i}_hindcast_success.json) )
$(A02_targets) : $(work_folder)/$(batch_key)/A02_prior/A02_%_hindcast_success.json : $(A01b_path)/A01b_ID_%.json $(analysisfolder)/A02b_WW3_hindcast_prior.py
					python $(analysisfolder)/A02b_WW3_hindcast_prior.py $* $(batch_key) $(ID_flag) > log/A02/$*.txt 2>&1

#/Users/Shared/Projects/2021_IceSAT2_tracks/data/scratch//SH_batch02/processed_ATL10_20190207001112_06190210_004_01.h5
.PHONY : A02 A03
A02 : $(A02_targets)
#A03 : $(A03_targets)


# A03_targets := $(foreach i, $(ID_track_list) ,$(addprefix $(A03_path)/processed_ATL10_, ${i}.h5 ) )
# $(A03_targets) : $(A03_path)/processed_ATL10_%.h5 : $(work_folder)B02_spectra_$(hemis)/B02_%_FFT.nc
# 					python $(track_downloader)/nsidc_icesat2_associated.py --user mhell@ucsd.edu --netrc ~/.netrc --product ATL03 --directory $(scratch_folder)/$(batch_key) -F ATL03_$(subst _003_,_005_,$(subst _004_,_005_,$*)).h5
# 					mv $(A03_path)/ATL10_$(subst _003_,_005_,$(subst _004_,_005_,$*)).h5 $(A03_path)/processed_ATL10_$(subst _003_,_005_,$(subst _004_,_005_,$*)).h5



### ---------- section B -------------
### derive beam list based on downloads
#A01c_success= $(shell ls $(A01c_download)/ATL03_*.h5)
#ID_track_list= $(foreach i, $(basename $(A01c_success)), $(subst $(A01c_download)/ATL03_,,$(i))  )

B01_targets := $(foreach i, $(ID_track_list) ,$(addprefix $(work_folder)/$(batch_key)/B01_regrid/, ${i}_B01_binned.h5 ) )
B02_targets := $(foreach i, $(ID_track_list) ,$(addprefix $(work_folder)/$(batch_key)/B02_spectra/, B02_${i}_FFT.nc ) )

.PHONY : B01 B02# B03
B01 : $(B01_targets)
B02 : $(B02_targets) #B02_collect

# !!!!!!!!! this does not work because the track ID chanegs from ALT10 to ALT03
$(B01_targets) : $(work_folder)/$(batch_key)/B01_regrid/%_B01_binned.h5 : $(A01b_path)/A01b_ID_%.json #$(analysisfolder)/B01_filter_regrid_segments.py
					python $(analysisfolder)/B01_filter_regrid_segments.py $* $(batch_key) $(ID_flag) > log/B01/$*.txt 2>&1

$(B02_targets) : $(work_folder)/$(batch_key)/B02_spectra/B02_%_FFT.nc : $(work_folder)/$(batch_key)/B01_regrid/%_B01_binned.h5 #$(analysisfolder)/B02_make_spectra.py
					sleep $${RANDOM:0:2}s
					python $(analysisfolder)/B02_make_spectra_gFT.py $* $(batch_key) $(ID_flag) > log/B02/$*.txt 2>&1

# B02_collect : #B01 B02
# 					${MKDIR_P} $(plot_folder)$(hemis)/$(batch_key)/B03/
# 					rm -fv $(plot_folder)$(hemis)/$(batch_key)/B03/*
# 					cp $(plot_folder)$(hemis)/$(batch_key)/*/B03_specs_*.png $(plot_folder)$(hemis)/$(batch_key)/B03/

#### delete bad tracks ###

# defeins targets for all bad tracks.
bad_beams_files := $(shell ls $(bad_tracks_folder)/$(batch_key)) #$(shell jq -r .[] $(config_folder)config.json)
bad_beams := $(basename $(bad_beams_files))

source_files := $(foreach i, $(bad_beams) ,$(addprefix $(scratch_folder)/$(batch_key), /ATL03_${i}.h5 ) )
bad_ID_files     := $(foreach i, $(bad_beams) ,$(addprefix $(A01b_path)/, A01b_ID_${i}.json) )x
#B01_files_02 := $(foreach i, $(bad_beams) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_corrected.h5 ) )
#B01_files_03 := $(foreach i, $(bad_beams) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_regridded.h5 ) )

plot_folders := $(foreach i, $(bad_beams) ,$(addprefix $(plot_folder)/$(hemis)/tracks/, ${i} ) )

rm_bad_beams: #rm_B01_files
					@echo "removing bad tracks in 5 sec"
					sleep 5
					rm -fv $(source_files)
					rm -fv -r $(plot_folders)
					rm -fr $(bad_ID_files)
					mv $(bad_tracks_folder)/$(batch_key)/*.json $(bad_tracks_folder)/$(batch_key)_deleted/

#rm -fv $(B01_targets) $(B01_files_02) $(B01_files_03)


# B03 overview plots
B02_path := $(work_folder)/$(batch_key)/B02_spectra/B02_
# find sources
B02_success := $(basename $(shell ls $(B02_path)*_gFT_x.nc) )
# B02_fail := $(basename $(shell ls $(B02_path)*_fail.json) )
B03_path := $(plot_folder)$(hemis)/$(batch_key)/
# write targets
B03_list := $(foreach i, $(B02_success) , $(subst $(B02_path),$(B03_path),$(i))  )
B03_targets := $(foreach i, $(B03_list) , $(subst _gFT_x,/B03_success.json,$(i) ) )


B03 : $(B03_targets) B03_mkdir B03_collect

reset_B03 :
					ls $(B03_path)*/B03_success.json
					rm -rf $(B03_path)*/B03_success.json
#.PHONY : B03_collect

B03_mkdir :
					${MKDIR_P} $(B03_collect_path)
					${MKDIR_P} $(B03_collect_path)coord_check

$(B03_targets) : $(B03_path)%/B03_success.json : $(work_folder)/$(batch_key)/B02_spectra/B02_%_FFT.nc
					python $(analysisfolder)/B03_plot_spectra_ov.py $* $(batch_key) $(ID_flag) > log/B03/$*.txt 2>&1

B03_collect_path := $(plot_folder)$(hemis)/$(batch_key)/B03_ov/
# search for all file and move, and rename them similar files
B03_success := $(shell ls $(B03_path)*/B03_success.json)
B03_collect_list := $(foreach i, $(B03_success) , $(subst $(B03_path),$(B03_collect_path),$(i) )  )
#B03_collect_targets := $(foreach i, $(B03_collect_list) , $(subst /B03_success.json,_B03_success.json,$(i)) )
B03_collect_targets := $(foreach i, $(B03_collect_list) , $(subst /B03_success.json,_B03_specs_L25000.png,$(i)) )
B03_collect : $(B03_collect_targets)

$(B03_collect_targets) : $(B03_collect_path)%_B03_specs_L25000.png : $(B03_path)%/B03_success.json
					cp $(B03_path)$*/B03_specs_L25000.0.png $(B03_collect_path)$*_B03_specs_L25000.png
					cp $(B03_path)$*/B03_specs_coord_check.png $(B03_collect_path)coord_check/$*_B03_specs_coord_check.png
#cp $(B03_path)$*/B03_spectra/B03_freq_reconst*.pdf $(B03_collect_path)$*_B03_freq_reconst.pdf



B04_targets := $(foreach i, $(B03_success) , $(subst B03_success,B04_success,$(i)) )

B04 : $(B04_targets) B04_collect

$(B04_targets) : $(B03_path)%/B04_success.json : $(B03_path)%/B03_success.json $(work_folder)/$(batch_key)/A02_prior/A02_%_hindcast_success.json
					python $(analysisfolder)/B04_angle.py $* $(batch_key) $(ID_flag) > log/B04/$*.txt 2>&1

B04_success := $(shell ls $(B03_path)*/B04_success.json)
B04_B05_collect_path := $(plot_folder)$(hemis)/$(batch_key)/B04_B05_angle/
B04_collect_list := $(foreach i, $(B04_success) , $(subst $(B03_path),$(B04_B05_collect_path),$(i) )  )
B04_collect_targets := $(foreach i, $(B04_collect_list) , $(subst /B04_success.json,_B04_marginal_distributions.pdf,$(i)) )

B04_B05_mkdir :
					${MKDIR_P} $(B04_B05_collect_path)

B04_collect : B04_B05_mkdir $(B04_collect_targets)
$(B04_collect_targets) : $(B04_B05_collect_path)%_B04_marginal_distributions.pdf : $(B03_path)%/B04_success.json
				cp $(B03_path)$*/B04_marginal_distributions.pdf $(B04_B05_collect_path)$*_B04_marginal_distributions.pdf

## ---------------- B5

B05_targets := $(foreach i, $(B04_success) , $(subst B04_success,B05_success,$(i)) )

B05 : $(B05_targets) B05_collect

$(B05_targets) : $(B03_path)%/B05_success.json : $(B03_path)%/B04_success.json
					python $(analysisfolder)/B05_define_angle.py $* $(batch_key) $(ID_flag) > log/B04/$*.txt 2>&1

B05_success := $(shell ls $(B03_path)*/B05_success.json)
B05_collect_list := $(foreach i, $(B05_success) , $(subst $(B03_path),$(B04_B05_collect_path),$(i) )  )
B05_collect_targets := $(foreach i, $(B05_collect_list) , $(subst /B05_success.json,_B05_dir_ov.pdf,$(i)) )

B05_collect : B04_B05_mkdir $(B05_collect_targets)
$(B05_collect_targets) : $(B04_B05_collect_path)%_B05_dir_ov.pdf : $(B03_path)%/B05_success.json
				cp $(B03_path)$*/B05_dir_ov.pdf $(B04_B05_collect_path)$*_B05_dir_ov.pdf

## ---------------- B6

B06_targets := $(foreach i, $(B05_success) , $(subst B05_success,B06_success,$(i)) )

B06 : $(B06_targets) B06_collect

$(B06_targets) : $(B03_path)%/B06_success.json : $(B03_path)%/B05_success.json
					python $(analysisfolder)/B06_correct_separate_var.py $* $(batch_key) $(ID_flag) > log/B04/$*.txt 2>&1

B06_success := $(shell ls $(B03_path)*/B06_success.json)
B06_collect_path := $(plot_folder)$(hemis)/$(batch_key)/B06_correction/
B06_collect_list := $(foreach i, $(B06_success) , $(subst $(B03_path),$(B06_collect_path),$(i) )  )
B06_collect_targets := $(foreach i, $(B06_collect_list) , $(subst /B06_success.json,_B06_atten_ov.png,$(i)) )

B06_mkdir :
					${MKDIR_P} $(B06_collect_path)

B06_collect : B06_mkdir $(B06_collect_targets)
$(B06_collect_targets) : $(B06_collect_path)%_B06_atten_ov.png : $(B03_path)%/B06_success.json
				cp $(B03_path)$*/B06_correction/$*_B06_atten_ov.png $(B06_collect_path)$*_B06_atten_ov.png


# ---------------- section C ------------


# ---------------- misc -----------------
# sync plots from batch key to gdrive folder.
sync_gdrive :
					rclone sync $(plot_folder)$(hemis)/$(batch_key)/ gdrive_brown:2021_ICESat2_tracks/plots/$(hemis)/$(batch_key)/
					#rclone sync $(plot_folder)../movie_temp_files/$(key) gdrive:Projects/2020_moist_two_layer/plots/movie_temp_files/$(key)/

# for printing variables, used for debugging
print-%  : ; @echo $* = $($*)


# dummy for sending mails when target is done.
## generated logfile at log/station.poc_log.txt
# ts := $(shell /bin/date "+%Y-%m-%d--%H-%M-%S")
#
# SUBJECT	=	"fig-dispoint - job $(ID)"
# EMAIL3	=	mhell@ucsd.edu
# mailer	=	/bin/mail -s $(SUBJECT) $(EMAIL3) < "$(logfile)"
#
# .PHONY : sendmail
# sendmail :
# 					@echo "---- Sendmail ----" >> $(logfile)
# 					$(mailer)
