MAKEFLAGS += -j12

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




analysisfolder		= /home/mhell/2021_ICESat2_tracks/analysis_db/
plotsfolder	     	= /home/mhell/2021_ICESat2_tracks/plots/

work_folder				= /work/mhell_work/2021_ICESat2_tracks/
scratch_folder    = /scratch/mhell/2021_ICESat2_tracks/
#movie_temp_folder	    = /home/mhell/2020_moist_two_layer/plots/movie_temp_files/
#targetfolder	     		= /home/mhell/2020_moist_two_layer/targets/
#codefolder			= /home/mhell/2020_moist_two_layer/code/

hemis 						= SH
batch							= batch01
#replace_flag			= False
beam 							= False

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

#### section B ####
B01_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_binned.h5 ) )
B02_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder)/B02_spectra_SH/, B02_${i}_FFT.nc ) )

.PHONY : B01 B02
B01 : $(B01_targets)
B02 : $(B02_targets)

$(B01_targets) : $(work_folder)/B01_regrid_$(hemis)/%_B01_binned.h5 : $(scratch_folder)/$(batch_key)/processed_ATL03_%.h5
					python $(analysisfolder)/B01_filter_regrid.py $* $(batch_key) $(test_flag) > log/B01/$*.txt 2>&1

$(B02_targets) : $(work_folder)/B02_spectra_$(hemis)/B02_%_FFT.nc : $(work_folder)/B01_regrid_$(hemis)/%_B01_binned.h5
					python $(analysisfolder)/B02_make_spectra.py $* $(batch_key) $(test_flag) > log/B02/$*.txt 2>&1


#
# sync_gdrive :
# 					rclone sync $(plotsfolder)$(key) gdrive:Projects/2020_moist_two_layer/plots/experiments/$(key)/
# 					#rclone sync $(plotsfolder)../movie_temp_files/$(key) gdrive:Projects/2020_moist_two_layer/plots/movie_temp_files/$(key)/


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
