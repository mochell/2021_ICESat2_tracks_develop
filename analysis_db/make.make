MAKEFLAGS += -j10

### example : show example how to use this file
.PHONY : example
example :
					@echo "to generate parameter collection "
					@echo "make -f run_model.make make_collection key=channel_N128_collection"
					@echo "then to run the models"
					@echo "make -f run_model.make run_models get_blocks key=channel_N128_collection1"
					@echo "Dry run, selecting a subset of year"
					@echo 'make -n -f master.make run_model.make run_models get_blocks key=channel_N128_collection1'
					@echo "force target to rebiuld"
					@echo "make -n -B -f master.make run_model.make run_models"


MKDIR_P = mkdir -p

analysisfolder		= /home/mhell/2021_ICESat2_tracks/analysis_db/
plotsfolder	     	= /home/mhell/2021_ICESat2_tracks/plots/

work_folder				= /work/mhell_work/2021_ICESat2_tracks/
scratch_folder    = /scratch/mhell/2021_ICESat2_tracks/
track_lists			  = /home/mhell/2021_ICESat2_tracks/track_lists/
track_downloader  = /home/mhell/2021_ICESat2_tracks/modules/read-ICESat-2/scripts/
#movie_temp_folder	    = /home/mhell/2020_moist_two_layer/plots/movie_temp_files/
#targetfolder	     		= /home/mhell/2020_moist_two_layer/targets/
#codefolder			= /home/mhell/2020_moist_two_layer/code/
bad_tracks_folder = $(work_folder)bad_tracks


### parameters
hemis = SH
batch = batch02#batch01 # define batch here. defined higher-level order structure.
#replace_flag			= False
beam 	= False

download=False # if not False provide a filename.json read this file and download the data. check that you are on the right batch!!

batch_key = $(hemis)_$(batch)
test_flag	= False
### ------

ifeq ($(beam), False)
	# load nc file list
	beam_list_json := $(shell ls $(scratch_folder)/$(batch_key)/) #$(shell jq -r .[] $(config_folder)config.json)

	beam_list_rar := $(basename $(beam_list_json))
	beam_list := $(foreach i, $(beam_list_rar) , $(subst processed_ATL03_,,$(i))  )

	#substr_list := $(foreach i, $(beam_list_rar) , $(subst _,'',$(i))  )
	#ALT := $(foreach i, $(beam_list_rar) , $(words $(subst _,'',$(i))  ) )
	# ALT := $(suffix $(beam_list_json) )
else
	beam_list= $(beam)
endif


# beam_list_rar := $(basename $(beam_list_json))
# beam_list := $(foreach i, $(beam_list_rar) , $(subst processed_ATL03_,,$(i))  )

.PHONY : make_folder
make_folder :
					${MKDIR_P} $(scratch_folder)/$(batch_key)
					${MKDIR_P} $(bad_tracks_folder)/$(batch_key)
					${MKDIR_P} $(bad_tracks_folder)/$(batch_key)_deleted


#### section A
# download data
ifeq ($(download), False)
A01_targets := $(foreach i, $(beam_list) ,$(addprefix $(scratch_folder)/$(batch_key)/, processed_ATL03_${i}.h5 ) )
$(A01_targets) : $(scratch_folder)/$(batch_key)/processed_ATL03_%.h5 :
					touch $(scratch_folder)/$(batch_key)/processed_ATL03_$*.h5

else

A01_track_names := $(shell jq -r .[] $(track_lists)$(download))
A01_tracks_rar := $(basename $(A01_track_names))
A01_tracks := $(foreach i, $(A01_tracks_rar) , $(subst ATL03_,,$(i))  )

A01_targets := $(foreach i, $(A01_tracks) ,$(addprefix $(scratch_folder)/$(batch_key)/, processed_ATL03_${i}.h5 ) )

$(A01_targets) : $(scratch_folder)/$(batch_key)/processed_ATL03_%.h5 :
					python $(track_downloader)/nsidc_icesat2_associated.py --user mhell@ucsd.edu --netrc ~/.netrc --product ATL03 --directory $(scratch_folder)/$(batch_key) -F ATL03_$*.h5
					mv $(scratch_folder)/$(batch_key)/ATL03_$*.h5 $(scratch_folder)/$(batch_key)/processed_ATL03_$*.h5
endif

.PHONY : A01
A01 : $(A01_targets)

#### section B ####
# Filter in and make spectra
ifeq ($(download), False)
beam_list := $(beam_list)
else
beam_list := $(A01_tracks)
endif

B01_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_binned.h5 ) )
B02_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder)/B02_spectra_SH/, B02_${i}_FFT.nc ) )

.PHONY : B01 B02
B01 : make_folder $(B01_targets)

B02 : $(B02_targets)

$(B01_targets) : $(work_folder)/B01_regrid_$(hemis)/%_B01_binned.h5 : $(scratch_folder)/$(batch_key)/processed_ATL03_%.h5
					python $(analysisfolder)/B01_filter_regrid.py $* $(batch_key) $(test_flag) > log/B01/$*.txt 2>&1

$(B02_targets) : $(work_folder)/B02_spectra_$(hemis)/B02_%_FFT.nc : $(work_folder)/B01_regrid_$(hemis)/%_B01_binned.h5
					python $(analysisfolder)/B02_make_spectra.py $* $(batch_key) $(test_flag) > log/B02/$*.txt 2>&1

B_collect : #B01 B02
					${MKDIR_P} $(plotsfolder)$(hemis)/$(batch_key)/B02/
					rm $(plotsfolder)$(hemis)/$(batch_key)/B02/*
					cp $(plotsfolder)$(hemis)/$(batch_key)/*/B_spectra/B02_specs_*.png $(plotsfolder)$(hemis)/$(batch_key)/B02/

#### delete bad tracks ###
bad_beams_files := $(shell ls $(bad_tracks_folder)/$(batch_key)) #$(shell jq -r .[] $(config_folder)config.json)
bad_beams := $(basename $(bad_beams_files))

source_files := $(foreach i, $(bad_beams) ,$(addprefix $(scratch_folder)/$(batch_key), /processed_ATL03_${i}.h5 ) )
#B01_files_01 := $(foreach i, $(bad_beams) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_binned.h5 ) )
B01_files_02 := $(foreach i, $(bad_beams) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_corrected.h5 ) )
B01_files_03 := $(foreach i, $(bad_beams) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_regridded.h5 ) )

plot_folders := $(foreach i, $(bad_beams) ,$(addprefix $(plotsfolder)/$(hemis)/tracks/, ${i} ) )

rm_bad_beams: #rm_B01_files
					@echo "removing bad tracks in 5 sec"
					sleep 5
					rm -fv $(source_files)
					rm -fv $(B01_targets) $(B01_files_02) $(B01_files_03)
					rm -fv -r $(plot_folders)
					mv $(bad_tracks_folder)/$(batch_key)/*.json $(bad_tracks_folder)/$(batch_key)_deleted/


# #beam_list := $(foreach i, $(beam_list_rar) , $(subst processed_ATL03_,,$(i))  )

sync_gdrive :
					rclone sync $(plotsfolder)$(hemis)/$(batch_key)/B02/ gdrive:Projects/2021_ICESat2_tracks/plots/$(hemis)/$(batch_key)/B02/
					#rclone sync $(plotsfolder)../movie_temp_files/$(key) gdrive:Projects/2020_moist_two_layer/plots/movie_temp_files/$(key)/


## generated logfile at log/station.poc_log.txt
ts := $(shell /bin/date "+%Y-%m-%d--%H-%M-%S")

SUBJECT	=	"fig-dispoint - job $(ID)"
EMAIL3	=	mhell@ucsd.edu
mailer	=	/bin/mail -s $(SUBJECT) $(EMAIL3) < "$(logfile)"

.PHONY : sendmail
sendmail :
					@echo "---- Sendmail ----" >> $(logfile)
					$(mailer)

print-%  : ; @echo $* = $($*)
