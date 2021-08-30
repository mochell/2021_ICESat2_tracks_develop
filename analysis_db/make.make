MAKEFLAGS += -j5

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
B1_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder)/B01_regrid_SH/, ${i}_B01_binned.h5 ) )

.PHONY : B01
B01 : $(B1_targets)

$(B1_targets) : $(work_folder)/B01_regrid_SH/%_B01_binned.h5 : $(scratch_folder)/$(batch_key)/processed_ATL03_%.h5
					python $(analysisfolder)/B01_filter_regrid.py $* $(batch_key) $(test_flag) > log/B01/$*.txt 2>&1




#collection_master:=$(codefolder)par_collections/$(key).json
#collection_targets:=$(addprefix $(codefolder)par_collections/$(key)/,$(beam_list_json) )
#collection_target:=$(par_collections)/$(key)
#run_targets	:= $(addprefix $(work_folder)/, ${beam_list})
#run_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i} ) )

#run_targets := $(foreach i, $(beam_list) ,$(addprefix $(scratch_folder),${i}/model_output.nc ) )
remove_targets := $(foreach i, $(beam_list) ,$(addprefix $(temp_folder_analysis_output), ${i}/a.nc ) )

run_targets := $(foreach i, $(beam_list) ,$(addprefix $(temp_folder_analysis_output), ${i}/model_output-invarients.nc ) )


run_test_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/model_output_basic.nc ) )
run_test_stokes_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/model_output_basic_stokes.nc ) )

formatted_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/model_output_formatted.nc ) )
spectra_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/power_spec.nc ) )

blocking_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/blocks_objects.nc ) )
blocking_time_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/blocks_time_stats.h5 ) )
LWA_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/LWA_blocks_objects.nc ) )
LWA_decor_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/LWA_measures.nc ) )

PV_budget_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/PV_budget.nc ) )
ENST_budget_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/pseudo_budget.nc ) )
rb_targets := $(foreach i, $(beam_list) ,$(addprefix $(work_folder),${i}/Rbeta_stats.nc ) )


merid_mean_targets := $(foreach i, $(beam_list) ,$(addprefix $(targetfolder)B01/,B01_${i}.json ) )
PV_grad_targets := $(foreach i, $(beam_list) ,$(addprefix $(targetfolder)B01/,B01_pv_grads_${i}.json ) )

B02_targets := $(foreach i, $(beam_list) ,$(addprefix $(targetfolder)B02/,B02_${i}.json ) )


C01_summary1_targets := $(foreach i, $(beam_list) ,$(addprefix $(targetfolder)C01/,C01_${i}.json ) )
C02_targets := $(foreach i, $(beam_list) ,$(addprefix $(targetfolder)C02/,C02_${i}.json ) )
C03_targets := $(foreach i, $(beam_list) ,$(addprefix $(targetfolder)C03/,C03_${i}.json ) )
C04_targets := $(foreach i, $(beam_list) ,$(addprefix $(targetfolder)C04/,C04_${i}.json ) )
C04_vid_file_targets := $(foreach i, $(beam_list) ,$(addprefix $(movie_temp_folder),${i}/ov_movie_f2_r10.mp4) )

.PHONY : run_model run_test_model
run_model : $(run_targets) $(formatted_targets)
run_model_only : $(run_targets)
run_test_model : $(run_test_stokes_targets) # $(run_test_targets)


get_blocks : $(blocking_targets) #$(blocking_time_targets)
get_budgets : $(PV_budget_targets) $(ENST_budget_targets)
get_spectra : $(spectra_targets)
get_stats : $(merid_mean_targets) $(blocking_targets) $(rb_targets) #$(blocking_time_targets)
get_rb : $(rb_targets) $(plotsfolder)/$(key)/$(key)_RB_binder.pdf
get_LWA : $(LWA_targets) $(LWA_decor_targets)

make_collection :
					python $(codefolder)/par_collection_creator.py $(key) $(replace_flag) $(cross_exploration)


$(run_targets) : $(temp_folder_analysis_output)%/model_output-invarients.nc : $(collection_folder)/%.json
					nice -n 5 python $(codefolder)/$(model_run_file) $* $(col_flag_pass) > code/logs/collections/$*.txt 2>&1
					mkdir -p $(temp_folder_analysis_output)$*
					cp $(temp_storage)model_output-$*_invarients.nc $(temp_folder_analysis_output)$*/model_output-invarients.nc
					sleep 120

#python $(codefolder)/moist_NL_QG_channel_optimised_parseeker_adaptive_time_budgets_nograds.py $* $(col_flag_pass) > code/logs/collections/$*.txt 2>&1

.PHONY : rerun_postprocessing reset
reset: reset_analysis remove_run
reset_analysis:
					rm -rf $(PV_grad_targets) $(merid_mean_targets) $(rb_targets)
					rm -rf $(C01_summary1_targets)
					rm -rf $(C02_targets)
					rm -rf $(C03_targets)
					rm -rf $(C04_targets)
					rm -rf $(C04_vid_file_targets)
					rm -rf $(movie_temp_folder)$(key)*

remove_run: $(remove_targets)

$(remove_targets) : $(temp_folder_analysis_output)%/a.nc :
					rm -rfv $(work_folder)$(key)*
					rm -rfv $(temp_folder_analysis_output)$(key)*
					rm -rfv $(temp_storage)model_output-$(key)*.nc
					rm -rfv $(plotsfolder)$(key)/*/
remove_collection:
					rm -rf $(collection_targets)
					rm -rf $(collection_folder)


rerun_targets := $(foreach i, $(beam_list) ,$(addprefix $(temp_folder_analysis_output), ${i}/rerun_analysis_log.json ) )
rerun_postprocessing : $(rerun_targets) $(formatted_targets)

$(rerun_targets) : $(temp_folder_analysis_output)%/rerun_analysis_log.json : $(temp_folder_analysis_output)%/model_output-invarients.nc
					python $(analysis_streamlined)/A00_rerun_postprocessing.py $* $(col_flag_pass) >> code/logs/collections/$*.txt 2>&1
					#cp $(temp_storage)model_output-%_invarients.nc $(temp_folder_analysis_output)$*/model_output-invarients.nc


$(formatted_targets) : $(work_folder)%/model_output_formatted.nc : $(temp_folder_analysis_output)%/model_output-invarients.nc #$(analysis_streamlined)/A01_reformat_data_and_plot_external.py
					sleep $$(( $$RANDOM%10))
					python $(analysis_streamlined)/A01_reformat_data_and_plot_external.py $* $(col_flag_pass) >> code/logs/collections/$*.txt 2>&1
					echo "formatting done"

# (formatted_targets) : $(work_folder)%/model_output_formatted.nc : $(scratch_folder)%/model_output.nc $(analysisfolder)/A00_reformat_data_and_plot.py
# python $(analysisfolder)/A00_reformat_data_and_plot.py $* $(col_flag_pass) >> code/logs/collections/$*.txt 2>&1


######
$(run_test_targets) : $(work_folder)%/model_output_basic.nc : $(collection_folder)/%.json
					python $(codefolder)/moist_NL_QG_channel_optimised_parseeker_adaptive_time_basics_hyper.py $* $(col_flag_pass) > code/logs/collections/$*_hyper.txt 2>&1

$(run_test_stokes_targets) : $(work_folder)%/model_output_basic_stokes.nc : $(collection_folder)/%.json
					python $(codefolder)/moist_NL_QG_channel_optimised_parseeker_adaptive_time_basics_stokes.py $* $(col_flag_pass) > code/logs/collections/$*_stokes.txt 2>&1


$(spectra_targets) : $(work_folder)%/power_spec.nc : $(work_folder)%/model_output_basic.nc $(analysisfolder)/A00_basic_spectra.py
					python $(analysisfolder)/A00_basic_spectra.py $* $(col_flag_pass) > code/logs/collections/A00/$*.txt 2>&1
######

$(blocking_targets) : $(work_folder)%/blocks_objects.nc : $(work_folder)%/model_output_formatted.nc $(analysis_streamlined)/A02_derive_blocking.py
					sleep $$(( $$RANDOM%100/10))
					python $(analysis_streamlined)/A02_derive_blocking.py $* $(col_flag_pass) > code/logs/collections/A02/$*.txt 2>&1


$(LWA_targets) : $(work_folder)%/LWA_blocks_objects.nc : $(work_folder)%/model_output_formatted.nc $(analysis_streamlined)/A02_derive_LWA.py
					python $(analysis_streamlined)/A02_derive_LWA.py $* $(col_flag_pass) > code/logs/collections/A02/$*.txt 2>&1


$(LWA_decor_targets) : $(work_folder)%/LWA_measures.nc : $(work_folder)%/model_output_formatted.nc $(analysis_streamlined)/S03_LWA_dist_decorr_pwlech.py
					python $(analysis_streamlined)/S03_LWA_dist_decorr_pwlech.py $* $(col_flag_pass) > code/logs/collections/A02/$*.txt 2>&1
# $(blocking_time_targets) : $(work_folder)%/blocks_time_stats.h5 : $(work_folder)%/model_output_formatted.nc $(analysisfolder)/A02_derive_blocking_tmean.py
# 					python $(analysisfolder)/A02_derive_blocking_tmean.py $* $(col_flag_pass) > code/logs/collections/A01/$*.txt 2>&1

$(PV_budget_targets) : $(work_folder)%/PV_budget.nc : $(work_folder)%/model_output_formatted.nc #$(analysisfolder)/A02_derive_PV_budget.py
					python $(analysisfolder)/A02_derive_PV_budget.py $* $(col_flag_pass) > code/logs/collections/A02/$*.txt 2>&1

$(ENST_budget_targets) : $(work_folder)%/pseudo_budget.nc : $(work_folder)%/model_output_formatted.nc $(analysis_streamlined)/A02_derive_enstrophy_budget.py
					python $(analysis_streamlined)/A02_derive_enstrophy_budget.py $* $(col_flag_pass) > code/logs/collections/A02/$*.txt 2>&1

$(rb_targets) : $(work_folder)%/Rbeta_stats.nc : $(work_folder)%/model_output_formatted.nc $(analysis_streamlined)/A02_derive_Rbeta.py
					python $(analysis_streamlined)/A02_derive_Rbeta.py $* $(col_flag_pass) > code/logs/collections/A02/$*.txt 2>&1


$(merid_mean_targets) : $(targetfolder)B01/B01_%.json : $(work_folder)%/blocks_objects.nc $(work_folder)%/pseudo_budget.nc $(analysis_streamlined)/B01_look_at_case_meridmean.py
					python $(analysis_streamlined)/B01_look_at_case_meridmean.py $* $(col_flag_pass) > code/logs/collections/B01/$*.txt 2>&1

$(PV_grad_targets) : $(targetfolder)B01/B01_pv_grads_%.json : $(work_folder)%/blocks_objects.nc $(work_folder)%/blocks_time_stats.h5 $(work_folder)%/pseudo_budget.nc $(analysisfolder)/B01_look_at_case_PV-grad.py
					python $(analysisfolder)/B01_look_at_case_PV-grad.py $* $(col_flag_pass) >> code/logs/collections/B01/$*.txt 2>&1

$(C01_summary1_targets) : $(targetfolder)C01/C01_%.json : $(work_folder)%/model_output_formatted.nc $(analysis_streamlined)/C01_summary1_Rb_composites.py $(work_folder)%/Rbeta_stats.nc $(work_folder)%/LWA_blocks_objects.nc $(work_folder)%/pseudo_budget.nc
					python $(analysis_streamlined)/C01_summary1_Rb_composites.py $* $(col_flag_pass) >> code/logs/collections/B01/$*.txt 2>&1

$(C02_targets) : $(targetfolder)C02/C02_%.json : $(work_folder)%/model_output_formatted.nc $(analysis_streamlined)/C02_blocking_index_hists.py $(work_folder)%/blocks_objects.nc $(work_folder)%/Rbeta_stats.nc $(work_folder)%/LWA_blocks_objects.nc $(analysis_streamlined)/C02_Rb_block_phase_hoffmoeller.py
					python $(analysis_streamlined)/C02_Rb_block_phase_hoffmoeller.py $* $(col_flag_pass) >> code/logs/collections/B01/$*.txt 2>&1
					python $(analysis_streamlined)/C02_blocking_index_hists.py $* $(col_flag_pass) >> code/logs/collections/B01/$*.txt 2>&1
					python $(analysis_streamlined)/C02_LWA_index_hists.py $* $(col_flag_pass) >> code/logs/collections/B01/$*.txt 2>&1 # this file makes the C02 target

$(C03_targets) : $(targetfolder)C03/C03_%.json : $(work_folder)%/model_output_formatted.nc $(analysisfolder)/C03_make_total_momentum_budget.py $(work_folder)%/Rbeta_stats.nc $(work_folder)%/LWA_blocks_objects.nc $(work_folder)%/pseudo_budget.nc
					python $(analysisfolder)/C03_make_total_momentum_budget.py $* $(col_flag_pass) >> code/logs/collections/B01/$*.txt 2>&1


#$(plotsfolder)/$(key)/$(key)_PV_grad_binder.pdf
make_binders: $(plotsfolder)/$(key)/$(key)_meridmean_binder.pdf  $(plotsfolder)/$(key)/$(key)_RB_binder.pdf $(plotsfolder)/$(key)/$(key)_summary_binder.pdf

make_summary: $(C01_summary1_targets) $(plotsfolder)/$(key)/$(key)_summary_binder.pdf
make_hist: $(C02_targets)
make_tot_mom : $(C03_targets) $(plotsfolder)/$(key)/$(key)_tot_mom_binder.pdf

$(plotsfolder)/$(key)/$(key)_meridmean_binder.pdf : $(merid_mean_targets)
					gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(plotsfolder)/$(key)/$(key)_meridmean_binder.pdf $(plotsfolder)/$(key)/meridmean/psi_ano_tmean_blocks_budget_*.pdf

$(plotsfolder)/$(key)/$(key)_PV_grad_binder.pdf : $(PV_grad_targets)
					gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(plotsfolder)/$(key)/$(key)_PV_grad_binder.pdf $(plotsfolder)/$(key)/PV_grad/pv_grads*.pdf

$(plotsfolder)/$(key)/$(key)_RB_binder.pdf : $(rb_targets)
					gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(plotsfolder)/$(key)/$(key)_RB_binder.pdf $(plotsfolder)/$(key)/Rbeta/rbeta_tmean_*.pdf

$(plotsfolder)/$(key)/$(key)_summary_binder.pdf : $(C01_summary1_targets) $(merid_mean_targets)
					gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(plotsfolder)/$(key)/$(key)_sumamry_binder.pdf $(plotsfolder)/$(key)/summaries/summary1_*.pdf
					#gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(plotsfolder)/$(key)/$(key)_homoeller_binder.pdf $(plotsfolder)/$(key)/summaries/RB_blocks_hoffmoeller_complete*.pdf
					gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(plotsfolder)/$(key)/$(key)_homoeller_binder_150d.pdf $(plotsfolder)/$(key)/summaries/RB_blocks_hoffmoeller_150_real_days*.pdf
					cp $(collection_master) $(plotsfolder)/$(key)/



$(plotsfolder)/$(key)/$(key)_tot_mom_binder.pdf : $(C03_targets)
					gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$(plotsfolder)/$(key)/$(key)_tot_mom_binder.pdf $(plotsfolder)/$(key)/C03_tot_mom/C03_tot_mombudget_hoffm_*.pdf $(plotsfolder)/$(key)/C03_tot_mom/C03_tot_mombudget_clim_*.pdf


#make_video: $(plotsfolder)/$(key)/$(key)_ov_movie_f2_r10.mp4
make_snapshots : $(C04_targets)
make_video : $(C04_vid_file_targets)
make_movie : make_video


$(C04_targets) : $(targetfolder)C04/C04_%.json : $(work_folder)%/LWA_blocks_objects.nc #$(work_folder)%/Rbeta_stats.nc
					python $(analysis_streamlined)/C04_map_movie_1colm.py $* $(col_flag_pass) >> code/logs/collections/B01/$*.txt 2>&1

ov_plots : $(B02_targets)
$(B02_targets) : $(targetfolder)B02/B02_%.json : $(work_folder)%/model_output_formatted.nc $(analysis_streamlined)/B02_plot_data.py #$(work_folder)%/Rbeta_stats.nc
					python $(analysis_streamlined)/B02_plot_data.py $* $(col_flag_pass) >> code/logs/collections/B01/$*.txt 2>&1

# $(plotsfolder)/$(key)/$(key)_ov_movie_f2_r10.mp4 : $(C04_targets)
# 					ffmpeg -y -framerate 2 -r 10 -f image2 -i $(plotsfolder)/$(key)/movies/temp_ov/movie_temp_%05d.png -f mp4 -q:v 0 -vcodec mpeg4 -r 10 $(plotsfolder)/$(key)/$(key)_ov_movie_f2_r10.mp4
# 					@echo "video done"

$(C04_vid_file_targets) : $(movie_temp_folder)%/ov_movie_f2_r10.mp4 :$(targetfolder)C04/C04_%.json
					ffmpeg -y -framerate 2 -r 10 -f image2 -i $(movie_temp_folder)$*/temp_files/movie_temp_%05d.png -f mp4 -q:v 0 -vcodec mpeg4 -r 10 $(movie_temp_folder)/$*_ov_movie_f2_r10.mp4
					@echo "video done"


# @if ($(collection_flag), True)
# 	# load json file
# 	beam_list_json := $(shell ls $(par_collections)/$(key)/) #$(shell jq -r .[] $(config_folder)config.json)
# 	beam_list := $(basename $(beam_list_json))
#
# else
# 	beam_list= $(key)
# endif
#
# @if [ "$(collection_flag)" = "True" ]; then\
# 	echo "collection_flag True";\
# 	python $(Bpath)/B01_fit_model_parallel.py $(key)
# else \
# 	echo "collection_flag False";\
# 	python $(Apath)/A01_load_NDBC_single.py $* $(config) >> log/A01/$*.nc.txt 2>&1;\
# fi

sync_gdrive :
					rclone sync $(plotsfolder)$(key) gdrive:Projects/2020_moist_two_layer/plots/experiments/$(key)/
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
