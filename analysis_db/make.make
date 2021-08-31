#MAKEFLAGS += -j12

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
#movie_temp_folder	    = /home/mhell/2020_moist_two_layer/plots/movie_temp_files/
#targetfolder	     		= /home/mhell/2020_moist_two_layer/targets/
#codefolder			= /home/mhell/2020_moist_two_layer/code/
bad_tracks_folder    = $(work_folder)bad_tracks

hemis 						= SH
batch							= batch01
#replace_flag			= False
beam 							= False

download 					= True

batch_key 				= $(hemis)_$(batch)
test_flag					= False


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

download_beam_list_json := $(shell ls $(scratch_folder)/$(batch_key)/) #$(shell jq -r .[] $(config_folder)config.json)

# beam_list_rar := $(basename $(beam_list_json))
# beam_list := $(foreach i, $(beam_list_rar) , $(subst processed_ATL03_,,$(i))  )

.PHONY : make_folder
make_folder :
					${MKDIR_P} $(bad_tracks_folder)/$(batch_key)
					${MKDIR_P} $(bad_tracks_folder)/$(batch_key)_deleted


#### section A ####
# download data

#A01_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_binned.h5 ) )

A01_prepare_data : A01_targets

ifeq ($(download_beams), True)

$(A01_targets) : $(scratch_folder)/$(batch_key)/processed_ATL03_%.h5
					python nsidc_icesat2_associated.py --user mhell@ucsd.edu --netrc ~/.netrc --product ATL03 ATL07-02_20200520190502_08440710_003_01.h5

else
$(A01_targets) : $(scratch_folder)/$(batch_key)/processed_ATL03_%.h5
					touch $(scratch_folder)/$(batch_key)/processed_ATL03_$*.h5

endif



#### section B ####
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
					${MKDIR_P} $(plotsfolder)$(hemis)/B02/
					rm $(plotsfolder)$(hemis)/B02/*
					cp $(plotsfolder)$(hemis)/tracks/*/B_spectra/B02_specs_*.png $(plotsfolder)$(hemis)/B02/

#### delete bad tracks ###
bad_beams_files := $(shell ls $(bad_tracks_folder)/$(batch_key)) #$(shell jq -r .[] $(config_folder)config.json)
bad_beams := $(basename $(bad_beams_files))

source_files := $(foreach i, $(bad_beams) ,$(addprefix $(scratch_folder)/$(batch_key), /processed_ATL03_${i}.h5 ) )
B01_files_01 := $(foreach i, $(bad_beams) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_binned.h5 ) )
B01_files_02 := $(foreach i, $(bad_beams) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_corrected.h5 ) )
B01_files_03 := $(foreach i, $(bad_beams) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_regridded.h5 ) )

plot_folders := $(foreach i, $(bad_beams) ,$(addprefix $(plotsfolder)/$(hemis)/tracks/, ${i} ) )

rm_bad_beams:
					@echo "removing bad tracks in 5 sec"
					sleep 5
					rm $(source_files)
					rm $(B01_files_01) $(B01_files_02) $(B01_files_03)
					rm -r $(plot_folders)
					mv $(bad_tracks_folder)/$(batch_key)/*.json $(bad_tracks_folder)/$(batch_key)_deleted/

# #beam_list := $(foreach i, $(beam_list_rar) , $(subst processed_ATL03_,,$(i))  )

sync_gdrive :
					rclone sync $(plotsfolder)$(hemis)/B02 gdrive:Projects/2021_ICESat2_tracks/plots/$(hemis)/
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
