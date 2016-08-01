function oned_model
%% ONED_MODEL simple ensemble data assimilation example.
%
%      There are no input arguments. Simply invoke by typing the name.
%
%      ONED_MODEL demonstrates the simplest possible case of ensemble data
%      assimilation. It is possible to explore assimilation algorithms,
%      ensemble sizes, model biases, etc. on-the-fly. The posterior
%      of the state is indicated by blue asterisks, the states evolve along
%      a trajectory indicated by the green lines to wind up at a prior state
%      for the assimilation - indicated by the green asterisks. After the
%      assimilation, the (posterior) state is indicated in blue and the
%      process is ready to repeat.
%
%      ONED_MODEL opens a gui control window that plots
%      the most recent prior, posterior, and observation, 
%      time sequences of the assimilation, the RMS error,
%      spread and kurtosis, and prior and posterior rank histograms.
%
%      The top button alternates between "Advance Model" and "Assimilate" to
%      single-step the model. The "Start Auto Run" button is useful to watch
%      the system evolve and generate estimates from many assimilation cycles.
%
%      Since this is a 'perfect model' experiment, we know the true state,
%      the amount of noise added to the observations, etc.; so it is possible to
%      calculate the error of the ensemble in addition to the spread. The
%      Truth is not (in general) the same as the observation!
%
%      This also introduces the concept of the 'rank histogram'. With N
%      ensemble members, there will always be N+1 'bins' that encompass the
%      state (to the left of the lowest ensemble member, to the right of the
%      highest, and all the ones in-between). The rank histogram tracks the
%      number of times the truth lies in each bin. With a good ensemble,
%      (i.e. a good assimilation system) the true state will be
%      indistinguishable from any of the ensemble members. This results
%      in a flat rank histogram, given enough samples.
%
% See also: gaussian_product, oned_ensemble, twod_ensemble, run_lorenz_63,
%           run_lorenz_96

%% DART software - Copyright 2004 - 2016 UCAR. This open source software is
% provided by UCAR, "as is", without charge, subject to all terms of use at
% http://www.image.ucar.edu/DAReS/DART/DART_download
%
% DART $Id$

help oned_model

atts = stylesheet; % get the default fonts and colors

%% -----------------------------------------------------------------------------

% Specify the figure size in pixels. After that, all positions are
% specified as fractions (units = Normalized). That way, the objects
% scale proportionally as the figure gets resized.

figXmin   =  50; % The horizontal position of the entire figure, in pixels
figYmin   = 250; % The vertical   position of the entire figure, in pixels
figWidth  = 900; % The width of the entire figure, in pixels
figHeight = 600; % The height of the entire figure, in pixels

handles.figure = figure('Position', [figXmin figYmin figWidth figHeight], ...
    'Units', 'Pixels', ...
    'Name', 'oned_model');

set(handles.figure,'Color', atts.background);

%% -----------------------------------------------------------------------------

%Position has units of normalized and therefore the components must be
%Fractions of figure. Also, the actual text  has FontUnits which are
%normalized, so the Font Size is a fraction of the text box/ edit box that
%contains it. Making the font size and box size normalized allows the text
%size and box size to change proportionately if the user resizes the
%figure window. This is true for all text/edit boxes and buttons.
% Define some variables to help specify locations/ratios on the figures.

scaled_fontsize  = 0.5;                                % set fontsize

handles.ui_button_advance_model = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.125 0.905 0.182 0.064], ...
    'String', 'Advance Model', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal', ...
    'Callback', @step_ahead);

handles.ui_button_start_auto_run = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.125 0.819 0.182 0.064], ...
    'String', 'Start Auto Run' , ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal', ...
    'Callback', @auto_run_Callback);

%% -------------------------t---------------------------------------------------
%  Set up a parent container so we can move the one container around instead of
%  trying to manipulate the positions of all the components.

handles.ui_InputPanel = uipanel('Units','Normalized', ...
    'Position', [0.774 0.378 0.210 0.370], ...
    'BorderType', 'none', ...
    'BackgroundColor', atts.background);

delta = 0.02;
tXmin = delta;
tXwid = 0.70;
tYwid = 0.23;

tYmin = 1.0 - tYwid - delta;

scaled_fontsize = 0.50;
handles.ui_text_ens_size = uicontrol(handles.ui_InputPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.010 0.750 0.690 0.180], ...
    'String', 'Ens. Size' , ...
    'BackgroundColor', atts.background, ...
    'HorizontalAlignment','right', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal');

eXmin = tXmin + tXwid + delta;

handles.ui_edit_ens_size = uicontrol(handles.ui_InputPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [eXmin tYmin+0.03 0.257 0.18], ...
    'String', '4', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal', ...
    'Callback', @ens_size_Callback);
handles.ens_size = str2double(get(handles.ui_edit_ens_size,'String'));

tYmin = tYmin - tYwid - delta;
handles.ui_text_model_bias = uicontrol(handles.ui_InputPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.010 0.488 0.690 0.180], ...
    'String', 'Model Bias' , ...
    'BackgroundColor', atts.background, ...
    'HorizontalAlignment','right', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal');

handles.ui_edit_model_bias = uicontrol(handles.ui_InputPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [eXmin tYmin+0.02 0.257 0.18], ...
    'String', '0.0', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal', ...
    'Callback', @model_bias_Callback);
handles.model_bias = str2double(get(handles.ui_edit_model_bias,'String'));

tYmin = tYmin - tYwid - delta;
handles.ui_text_inflation = uicontrol(handles.ui_InputPanel, ...
    'Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.010 0.240 0.690 0.180], ...
    'String', 'Inflation' , ...
    'BackgroundColor', atts.background, ...
    'HorizontalAlignment','right', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal');

handles.ui_edit_inflation = uicontrol(handles.ui_InputPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [eXmin tYmin+0.02 0.257 0.18], ...
    'String', '1.0', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal', ...
    'Callback', @inflation_Callback);
handles.inflation = str2double(get(handles.ui_edit_inflation,'String'));

tYmin = tYmin - tYwid - delta;
handles.ui_text_nonlin_a = uicontrol(handles.ui_InputPanel, ...
    'Style', 'text',...
    'Units', 'Normalized',...
    'Position', [0.010 0.007 0.690 0.180],...
    'String', 'Nonlin. a' ,...
    'BackgroundColor', atts.background,...
    'HorizontalAlignment','right', ...
    'FontName', atts.fontname,...
    'FontUnits', 'Normalized',...
    'FontSize', scaled_fontsize,...
    'FontWeight', 'normal');

handles.ui_edit_nonlin_a = uicontrol(handles.ui_InputPanel, ...
    'Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [eXmin tYmin+0.03 0.257 0.18], ...
    'String', '0.0', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal', ...
    'Callback', @nonlin_a_Callback);
handles.nonlin_a = str2double(get(handles.ui_edit_nonlin_a,'String'));

%% -----------------------------------------------------------------------------

handles.ui_radio_button_group = uibuttongroup('BackgroundColor', atts.background, ...
    'BorderType', 'none', ...
    'Position',[0.834 0.142 0.130 0.240]);

rXmin    = 0.033;  % x location of radio buttons
rdeltaX  = 0.830;  % width of radio buttons
rHeight  = 0.275;  % height of radio buttons
dy       = (1.0 - 3*rHeight)/2.0;

rypos = 1.0 - rHeight - dy;
handles.ui_radio_button_eakf = uicontrol(handles.ui_radio_button_group, ...
    'Style', 'radio button', ...
    'Units', 'Normalized', ...
    'Position', [rXmin, rypos, rdeltaX , rHeight], ...
    'String', 'EAKF', ...
    'BackgroundColor', atts.background, ...
    'Foreground', 'k', ...
    'FontUnits', 'Normalized', ...
    'FontName', atts.fontname, ...
    'FontWeight','Bold', ...
    'FontSize', scaled_fontsize, ...
    'HandleVisibility', 'On');

rypos = rypos - rHeight;
handles.ui_radio_button_enkf = uicontrol(handles.ui_radio_button_group, ...
    'Style', 'radio button', ...
    'Units', 'Normalized', ...
    'Position', [rXmin, rypos, rdeltaX , rHeight], ...
    'String', 'EnKF' , ...
    'BackgroundColor', atts.background, ...
    'Foreground', 'k', ...
    'FontUnits', 'Normalized', ...
    'FontName', atts.fontname, ...
    'FontWeight','Bold', ...
    'FontSize', scaled_fontsize, ...
    'HandleVisibility', 'On');

rypos = rypos - rHeight;
handles.ui_radio_button_rhf = uicontrol(handles.ui_radio_button_group, ...
    'Style', 'radio button', ...
    'Units', 'Normalized', ...
    'Position', [rXmin, rypos, rdeltaX , rHeight], ...
    'String', 'RHF', ...
    'BackgroundColor', atts.background, ...
    'Foreground', 'k', ...
    'FontUnits', 'Normalized', ...
    'FontName', atts.fontname, ...
    'FontWeight','Bold', ...
    'FontSize', scaled_fontsize, ...
    'HandleVisibility', 'On');

%% -----------------------------------------------------------------------------

handles.ui_button_reset = uicontrol('Style', 'pushbutton', ...
    'Units', 'Normalized', ...
    'Position', [0.836 0.060 0.095 0.063], ...
    'String', 'Reset', ...
    'BackgroundColor', 'White', ...
    'FontName', atts.fontname, ...
    'FontUnits', 'Normalized', ...
    'FontSize', scaled_fontsize, ...
    'FontWeight', 'normal', ...
    'Callback', @reset_button_Callback);

%% -----------------------------------------------------------------------------
%  These appear to be error messages that can be turned on or off.
%  They start turned off.

handles.ui_text_inf_err_print = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.0600 0.6000 0.3000 0.0800], ...
    'String', 'ERROR: Inflation value must be between 1 and 5.', ...
    'ForegroundColor', atts.red, ...
    'BackgroundColor', 'White', ...
    'FontSize', atts.fontsize, ...
    'FontName', atts.fontname, ...
    'FontWeight', 'Bold', ...
    'Visible','Off');

handles.ui_text_model_bias_err_print = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.0600 0.6000 0.3000 0.0800], ...
    'String', 'ERROR: Model Bias value must be greater or equal to 0.', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.red, ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize, ...
    'FontWeight', 'Bold', ...
    'Visible','Off');

handles.ui_text_nonlin_err_print = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.0600 0.6000 0.3000 0.0800], ...
    'String', 'ERROR: Nonlin a must be under or equal to 0.', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.red, ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize, ...
    'FontWeight', 'Bold', ...
    'Visible','Off');

handles.ui_text_ens_size_err_print = uicontrol('Style', 'text', ...
    'Units', 'Normalized', ...
    'Position', [0.0600 0.6000 0.3000 0.0800], ...
    'String', 'ERROR: Ens. Size value must be greater or equal to 2.', ...
    'BackgroundColor', 'White', ...
    'ForegroundColor', atts.red, ...
    'FontName', atts.fontname, ...
    'FontSize', atts.fontsize, ...
    'FontWeight', 'Bold', ...
    'Visible','Off');

%% -----------------------------------------------------------------------------
%  Set the default values for everyone

reset_button_Callback()

%% -----------------------------------------------------------------------------

    function ens_size_Callback(~, ~)
        
        new_ens_size = str2double(get(handles.ui_edit_ens_size, 'String'));
        
        % Get a new ensemble size if not valid value
        if( ~ isfinite(new_ens_size) || (new_ens_size < 2) )
            
            fprintf('ERROR: Ens. Size value must be greater or equal to 2.\n')
            fprintf('ERROR: Ens. Size value must be greater or equal to 2.\n')
            
            % After this, only this edit box will work
            turn_off_controls;
            
            set(handles.ui_edit_ens_size, 'Enable', 'On',  ...
                                          'String', '?', ...
                                          'BackgroundColor', 'r');
            set(handles.ui_text_ens_size_err_print, 'Visible', 'On')
            
            return
        end
        
        turn_on_controls;
        
        set(handles.ui_edit_ens_size, 'Enable', 'On', 'BackgroundColor', 'White');
        set(handles.ui_text_ens_size_err_print, 'Visible', 'Off')
        
        % Generate a new ensemble by truncating old ensemble OR adding new
        if(new_ens_size == handles.ens_size)
            
            return
            
        elseif(new_ens_size < handles.ens_size)
            
            % Get rid of extra ensemble members, recompute mean, spread and kurtosis
            handles.ens      = handles.ens(1:new_ens_size);
            handles.ens_size = new_ens_size;
            
        else
            % Add new ensemble members drawn from present distribution
            handles.ens(handles.ens_size + 1 : new_ens_size) = ...
                randn([1 new_ens_size - handles.ens_size]);
            handles.ens_size = new_ens_size;
            
        end
        
        % Update moments
        handles.error    = calculate_rmse(handles.ens, 0.0);
        handles.spread   = std(handles.ens);
        handles.kurtosis = kurt(handles.ens);
        
        % If you change the ensemble size, you also have to reset the
        % histograms.
        handles.prior_rank = zeros(1, handles.ens_size + 1);
        handles.post_rank  = zeros(1, handles.ens_size + 1);
        
        set_prior_histogram();
        set_posterior_histogram();
        
    end

%% -----------------------------------------------------------------------------

    function model_bias_Callback(~, ~)
        
        % Check to make sure the input is a valid number
        model_bias_value = str2double(get(handles.ui_edit_model_bias, 'String'));
        
        if(isfinite(model_bias_value) && (model_bias_value >= 0))

            % If valid, update the value of the model bias. 
            handles.model_bias = model_bias_value;
            turn_on_controls;
            set(handles.ui_text_model_bias_err_print,'Visible','Off')
            set(handles.ui_edit_model_bias, 'Enable', 'On', ...
                                            'BackgroundColor', 'White');
            
        else
            % If not valid, force user to try again.
            % After this, only this edit box will work
            turn_off_controls;
            set(handles.ui_edit_model_bias, 'String', '?', ...
                                            'Enable', 'On', ...
                                            'BackgroundColor', 'r');
            set(handles.ui_text_model_bias_err_print,'Visible','On')
            
            fprintf('ERROR: Model Bias value must be greater or equal to 0.\n')
            fprintf('ERROR: Model Bias value must be greater or equal to 0.\n')
            
            return
        end
    end

%% -----------------------------------------------------------------------------

    function inflation_Callback(~, ~)
        
        % Get the value of the inflation
        inflation_value = str2double(get(handles.ui_edit_inflation, 'String'));
        
        if(isfinite(inflation_value) && (inflation_value >= 1) && (inflation_value <= 5))
            
            handles.inflation = inflation_value;
            
            turn_on_controls;
            
            set(handles.ui_edit_inflation, 'Enable', 'On', 'BackgroundColor', 'White');
            set(handles.ui_text_inf_err_print,'Visible','Off')
            
        else
            
            fprintf('ERROR: Inflation value must be between 1 and 5.\n')
            fprintf('ERROR: Inflation value must be between 1 and 5.\n')
            
            % After this, only this edit box will work
            turn_off_controls;
            
            set(handles.ui_edit_inflation, 'Enable', 'On', ...
                                           'String', '?', ...
                                           'BackgroundColor', 'r');
            set(handles.ui_text_inf_err_print,'Visible','On')
            
            return
            
        end
        
    end

%% -----------------------------------------------------------------------------

    function nonlin_a_Callback(~, ~)
        
        % Get the value of the model nonlinearity parameter 'alpha'
        
        nonlin_value = str2double(get(handles.ui_edit_nonlin_a, 'String'));
        
        if(isfinite(nonlin_value) && (nonlin_value <= 0))
            
            handles.alpha = nonlin_value;
            turn_on_controls;
            
            set(handles.ui_edit_nonlin_a, 'Enable', 'On', 'BackgroundColor', 'White');
            set(handles.ui_text_nonlin_err_print, 'Visible', 'Off')
            
        else  % ERROR STATE, force them to fix before moving on
            
            % After this, only this edit box will work
            turn_off_controls;
            
            fprintf('ERROR: Nonlin a must be under or equal to 0.\n')
            fprintf('ERROR: Nonlin a must be under or equal to 0.\n')
            
            set(handles.ui_edit_nonlin_a, 'Enable', 'On', ...
                                          'String', '?', ...
                                          'BackgroundColor', 'r');
            set(handles.ui_text_nonlin_err_print, 'Visible', 'On')
            
            return
            
        end
        
    end

%% -----------------------------------------------------------------------------


    function reset_button_Callback(~, ~)

        initialize_data();
        reset_graphics();
        
    end

%% -----------------------------------------------------------------------------


    function initialize_data(~, ~)
        
        % Reset all the figures and the data structures
        % Keep the current filter type, ensemble size and obs characteristics
        % Reset the time to 1 and be ready to advance
        
        % set random number seed to same value to generate known sequences
        % rng('default') is the Mersenne Twister with seed 0
        rng(0,'twister')
        
        % Set up global storage with initial values
        handles.ens_size         = 4;
        handles.ens              = randn(1, handles.ens_size);
        handles.model_bias       = 0.0;
        handles.inflation        = 1.0;
        
        handles.time_step        = 1;
        handles.ready_to_advance = true;
        handles.alpha            = 0.0;   % aka nonlin a
        handles.obs_error_sd     = 1;
        handles.observation      = 0;
        
        %  Compute the initial error (truth is 0.0) and spread (standard deviation)
        handles.error    = calculate_rmse(handles.ens, 0.0);
        handles.spread   = std(handles.ens);
        handles.kurtosis = kurt(handles.ens);
        
        % An array to keep track of rank histograms
        handles.prior_rank = zeros(1, handles.ens_size + 1);
        handles.post_rank  = zeros(1, handles.ens_size + 1);
        
        % Set the ui values/strings to starting values.
        set(handles.ui_button_advance_model, 'String', 'Advance Model');
        
        set(handles.ui_edit_ens_size,   'Value',               handles.ens_size);
        set(handles.ui_edit_ens_size,   'String', sprintf('%d',handles.ens_size));
        
        set(handles.ui_edit_model_bias, 'Value',                 handles.model_bias);
        set(handles.ui_edit_model_bias, 'String', sprintf('%.1f',handles.model_bias));
        
        set(handles.ui_edit_inflation,  'Value',                 handles.inflation);
        set(handles.ui_edit_inflation,  'String', sprintf('%.1f',handles.inflation));
        
        set(handles.ui_edit_nonlin_a,   'Value',                 handles.alpha);
        set(handles.ui_edit_nonlin_a,   'String', sprintf('%.1f',handles.alpha));
        
    end

%% -----------------------------------------------------------------------------

    function reset_graphics(~, ~)
        
        set_main_axes();
        set_error_spread_evolution();
        set_kurtosis_evolution()
        set_prior_histogram();
        set_posterior_histogram();
        set_state_evolution();
        
    end

%% -----------------------------------------------------------------------------

    function auto_run_Callback(~, ~)
        
        % Turn off all the other controls to avoid a mess
        turn_off_controls;
        
        set(handles.ui_button_start_auto_run, 'Enable', 'On');
        
        if(strcmp(get(handles.ui_button_start_auto_run, 'String'), 'Pause Auto Run'))
            
            % Being told to stop; switch to not running status
            set(handles.ui_button_start_auto_run, 'String', 'Start Auto Run');
            
        else
            % Being told to start run
            % Change the button to 'Pause Auto Run')
            set(handles.ui_button_start_auto_run, 'String', 'Pause Auto Run');
            
            % Loop through advance and assimilate steps until stopped
            while(true)
                
                % Check to see if stop has been pushed
                status_string = get(handles.ui_button_start_auto_run, 'String');
                
                if(strcmp(status_string, 'Start Auto Run'))
                    
                    turn_on_controls;
                    
                    return
                    
                end
                
                % Do the next advance or assimilation step
                step_ahead;
                pause(0.2)
                
            end
            
        end
        
        % Turn all the other controls back on
        turn_on_controls;
        
    end

%% -----------------------------------------------------------------------------

    function step_ahead(~, ~)
        
        % Start out working in ensemble time series plot
        axes(handles.h_state_evolution);
        
        % If this is an advance, get and plot new model state, advance time
        if(handles.ready_to_advance);
            
            % Set to do an assimilation next time
            handles.ready_to_advance = false;
            
            % Advance the model and then inflate
            ens_new = advance_oned(handles.ens, handles.alpha, handles.model_bias);
            ens_new_mean = mean(ens_new);
            ens_new = (ens_new - ens_new_mean) * sqrt(handles.inflation) + ens_new_mean;
            
            % plot the model evolution
            handles.time_step = handles.time_step + 1;
            h_evolution.prior = plot(handles.time_step - 0.1, ens_new, '*', ...
                'MarkerSize', 6, 'Color', atts.green);
            
            for i = 1:handles.ens_size
                plot([handles.time_step - 1 + 0.1, handles.time_step - 0.1], ...
                    [handles.ens(i), ens_new(i)], 'Color', atts.green);
            end
            
            %% Plot the segment for the prior error and spread
            % Want the lower y limit to stay 0 for error spread
            axes(handles.h_err_spread_evolution);
            
            prior_error = calculate_rmse(ens_new, 0.0);
            prior_spread = std(ens_new);
            
            h_e = line([handles.time_step - 1 + 0.1, handles.time_step - 0.1], ...
                [handles.error, prior_error]);
            set(h_e, 'Color', 'b', 'LineWidth', 2.0);
            
            h_s = line([handles.time_step - 1 + 0.1, handles.time_step - 0.1], ...
                [handles.spread, prior_spread]);
            set(h_s, 'Color', 'r', 'LineWidth', 2.0);
            
            handles.error  = prior_error;
            handles.spread = prior_spread;
            
            legend([h_e h_s], 'Error', 'Spread', 'Location', 'NorthEast');
            set(legend,'FontName', atts.fontname, 'FontSize', atts.fontsize);
            legend boxon
            
            axlims    = axis;
            axlims(3) = 0.0;
            axis(axlims)
            
            %% Plot the segment for the prior kurtosis
            % Want the lower y limit to stay 0 for kurtosis
            
            axes(handles.h_kurtosis_evolution);
            
            prior_kurtosis = kurt(ens_new);
            
            plot([handles.time_step - 1 + 0.1, handles.time_step - 0.1], ...
                [handles.kurtosis, prior_kurtosis], 'Color', 'r','LineWidth',2);
            
            handles.kurtosis = prior_kurtosis;

            axlims = axis;
            axlims(3) = 0.0;
            axis(axlims)
            
            % Update the prior rank histogram figure
            axes(handles.h_prior_rank_histogram);
            
            ens_rank = get_ens_rank(ens_new, 0);
            
            % Plot the latest rank entry as a different color
            temp_rank(:, 1)        = handles.prior_rank(1:handles.ens_size + 1);
            temp_rank(:, 2)        = 0;
            temp_rank(ens_rank, 2) = 1;
            
            bar(temp_rank,'stacked');
            
            % Plot the figure window for this update
            axes(handles.axes);
            hold off;
            
            % Find the limits of the plot
            % The height of the obs likelihood controls the vertical axis
            y_max = 1 / (sqrt(2 * pi) * handles.obs_error_sd);
            xmin = -10;
            xmax =  10;
            
            % Put on a black axis line using data limits

            plot([xmin xmax], [0, 0], 'k', 'Linewidth', 2);
            hold on;
            ens_axis = [xmin xmax -0.2 y_max + 0.02];
            
            % Set the axis for the figure window
            axis(ens_axis);
            
            % Turn off the negative labels
            set(gca, 'YTick', [0 0.1 0.2 0.3 0.4]);
            grid on;

            % Plot the prior ensemble members in green
            % Plotting ticks instead of asterisks makes bins clearer
            tick_half = 0.015;
            
            for n_tick = 1:handles.ens_size
                hg_prior = line([ens_new(n_tick), ens_new(n_tick)], ...
                    [-tick_half, tick_half]);
                set(hg_prior, 'Color', atts.green, 'LineWidth', 2);
            end
            
            % Plot the truth (at 0) as a tick
            hg_truth = line([0 0], [-0.02 0.02]);
            set(hg_truth, 'Color', 'k', 'LineWidth', 2);
            
            % Put in some information about the bin, x position tricky
            base_x = max(0, ens_axis(3));
            text_width = (ens_axis(4) - ens_axis(3)) / 3;
            
            if((base_x + text_width) > ens_axis(4))
                base_x = ens_axis(4) - text_width;
            end
            
            text(base_x, -0.1, ['Truth in Prior Bin ', num2str(ens_rank)], ...
                'FontSize', 14, 'FontWeight', 'Bold','FontName', atts.fontname);
            
            % Draw a line from the label string to the truth
            h = line([base_x + text_width / 8, 0], [-0.08, -0.03]);
            set(h, 'Color', 'k');
            
            % Label this plot
            xlabel('State'               ,'FontName', atts.fontname, 'FontSize', atts.fontsize);
            title('Latest Ensemble Prior','FontName', atts.fontname, 'FontSize', atts.fontsize);
            
            legend([hg_prior hg_truth],'Prior','Truth');
            set(legend                   ,'FontName', atts.fontname, 'FontSize', atts.fontsize);
            
            % Update the permanent storage of the rank values
            handles.prior_rank(ens_rank) = handles.prior_rank(ens_rank) + 1;
            
            % Set the pushbutton to say Assimilate Obs
            set(handles.ui_button_advance_model, 'String', 'Assimilate Obs');
            
            % Update the global storage of the ensemble
            handles.ens = ens_new;
            
        else
            % Ready to do an assimilation
            % Next step should be an advance
            handles.ready_to_advance = true;
            
            % Generate the observation as a draw Normal(0, 1)
            obs_error_sd = handles.obs_error_sd;
            observation = obs_error_sd * randn(1);
            
            % Plot the observation
            plot(handles.time_step, observation, 'r*', 'MarkerSize', 10);
            
            % Set the pushbutton to say Advance Model
            set(handles.ui_button_advance_model, 'String', 'Advance Model');
            
            % Adjust the horizontal range of the plot windows as needed
            if( mod(handles.time_step, 5) == 0)
                
                axes(handles.h_state_evolution);
                axlims    = axis;
                axlims(1) = handles.time_step - 4;
                axlims(2) = handles.time_step + 6;
                axis(axlims)
                
                % Want the lower y limit to stay 0 for error spread
                
                axes(handles.h_err_spread_evolution);
                axlims    = axis;
                axlims(1) = handles.time_step - 4;
                axlims(2) = handles.time_step + 6;
                axlims(3) = 0.0;
                axis(axlims)
                
                % Want the lower y limit to stay 0 for kurtosis
                
                axes(handles.h_kurtosis_evolution);
                axlims    = axis;
                axlims(1) = handles.time_step - 4;
                axlims(2) = handles.time_step + 6;
                axlims(3) = 0.0;
                axis(axlims)
                
            end
            
            % Do the assimilation
            ens = handles.ens;
            obs_error_sd = handles.obs_error_sd;
            
            % Figure out which filter option is currently selected
            val = get(handles.ui_radio_button_group,'SelectedObject');
            filter_type = get(val,'String');
            
            switch filter_type
                
                case 'EAKF'
                    [obs_increments, ~] = ...
                        obs_increment_eakf(ens, observation, obs_error_sd^2);
                case 'EnKF'
                    [obs_increments, ~] = ...
                        obs_increment_enkf(ens, observation, obs_error_sd^2);
                case 'RHF'
                    [obs_increments, ~] = ...
                        obs_increment_rhf(ens, observation, obs_error_sd^2);
            end
            
            new_ens = ens + obs_increments;
            axes(handles.h_state_evolution);
            plot(handles.time_step + 0.1, new_ens, 'b*', 'MarkerSize', 6);
            handles.ens = new_ens;
            
            % Update the rank data
            axes(handles.h_post_rank_histogram);
            ens_rank = get_ens_rank(handles.ens, 0);
            
            % Plot the latest rank entry as a different color
            temp_rank(:, 1)        = handles.post_rank(1:handles.ens_size + 1);
            temp_rank(:, 2)        = 0;
            temp_rank(ens_rank, 2) = 1;
            
            bar(temp_rank, 'stacked');
            
            % Update the permanent storage of the rank values
            handles.post_rank(ens_rank) = handles.post_rank(ens_rank) + 1;
            
            %% Plot the segment for the updated error
            
            axes(handles.h_err_spread_evolution);
            
            post_error = calculate_rmse(new_ens, 0.0);
            
            h = plot([handles.time_step - 0.1, handles.time_step + 0.1], ...
                [handles.error, post_error]);
            set(h,'Color', 'b', 'LineWidth', 2.0);
            
            handles.error = post_error;
            
            % Plot the segment for the updated spread
            post_spread = std(new_ens);
            
            axes(handles.h_err_spread_evolution);
            
            h = plot([handles.time_step - 0.1, handles.time_step + 0.1], ...
                [handles.spread, post_spread]);
            set(h, 'Color', 'r', 'LineWidth', 2.0);
            
            % Want the lower y limit to stay 0 for error spread
            axlims    = axis;
            axlims(3) = 0.0;
            axis(axlims)
            
            handles.spread = post_spread;
            
            %% Plot the segment for the updated kurtosis
            % Want the lower y limit to stay 0 for kurtosis
            
            axes(handles.h_kurtosis_evolution);
            
            post_kurtosis = kurt(new_ens);
            
            h = plot([handles.time_step - 0.1, handles.time_step + 0.1], ...
                [handles.kurtosis, post_kurtosis]);
            set(h, 'Color', 'r', 'LineWidth', 2.0);
            
            axlims    = axis;
            axlims(3) = 0.0;
            axis(axlims)
            
            handles.kurtosis= post_kurtosis;
            
            %% Plot the figure for this update
            axes(handles.axes);
            hold off
            
            % Find the limits of the plot
            % The height of the obs likelihood controls the vertical axis
            % Plot the observation likelihood
            [hg_like, ~, ylims] = plot_gaussian(observation, obs_error_sd, 1.0);
            set(hg_like, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
            hold on;
            
            xmin = -10;
            xmax =  10;
            ens_axis = [xmin xmax -0.2 ylims(2)+0.02];
            axis(ens_axis);
            
            % Put on a black axis line using data limits
            plot([xmin xmax], [0, 0], 'k', 'Linewidth', 2);

            % Plot the prior ensemble members in green
            % Plotting ticks instead of asterisks makes bins clearer
            tick_half = 0.015;
            
            for n_tick = 1:handles.ens_size
                hg_prior = plot([ens(n_tick), ens(n_tick)], ...
                    [-tick_half, tick_half], 'Color', atts.green, ...
                    'LineWidth', 2);
            end
            
            % Plot the posterior ensemble members in blue
            for n_tick = 1:handles.ens_size
                hg_post = plot([new_ens(n_tick), new_ens(n_tick)], ...
                    [-0.1 - tick_half, -0.1 + tick_half], 'Color', 'b', ...
                    'LineWidth', 2);
            end
                        
            % Plot the observation (at 0) as an asterisk
            plot(observation, 0, 'r*', 'MarkerSize', 14, 'LineWidth', 2.0);
            
            % Plot the truth (at -0.1 and 0) as a tick
            tick_half = 0.02;
            plot([0 0], [-0.1 - tick_half, -0.1 + tick_half], ...
                'k', 'LineWidth', 2.0);
            
            plot([0 0], [- tick_half, tick_half], ...
                'k', 'LineWidth', 2.0);
            
            % Put in some information about the bin, x position tricky
            base_x = max(0, ens_axis(1));
            text_width = (ens_axis(2) - ens_axis(1)) / 3;
            
            if((base_x + text_width) > ens_axis(2))
                base_x = ens_axis(2) - text_width;
            end
            
            text(base_x, -0.18, ['Truth in Posterior Bin ', num2str(ens_rank)], ...
                'FontSize', 14, 'FontWeight', 'Bold','FontName', atts.fontname);
            
            % Draw a line from the label string to the truth
            plot([base_x + text_width/8, 0], [-0.16, -0.13], 'k');
            
            % Fix up the final axis
            def_axis = axis;
            def_axis(1:2) = ens_axis(1:2);
            def_axis(3) = -0.2;
            def_axis(4) = 1 / (sqrt(2 * pi) * handles.obs_error_sd) + 0.02;
            axis(def_axis);
            
            % Turn off the negative labels
            set(gca, 'YTick', [0 0.1 0.2 0.3 0.4]);
            grid on;
            
            % Plot an additional axis
            plot(ens_axis(1:2), [-0.1 -0.1], 'k', 'LineWidth', 2);
            
            % Label this plot
            xlabel('State','FontName', atts.fontname,'FontSize', atts.fontsize);
            title('Latest Ensemble Prior, Likelihood, Posterior', ...
                'FontName', atts.fontname, 'FontSize', atts.fontsize,'FontWeight', 'Bold');
            
            % Put on legend
            legend([hg_prior hg_post hg_like], 'Prior', 'Posterior', 'Likelihood', 'Location','NorthEast');
            set(legend,'FontName', atts.fontname, 'FontSize', atts.fontsize);
        end
        
    end

%% -----------------------------------------------------------------------------

    function turn_off_controls()
        
        % Turn off all the other controls to avoid a mess
        set(handles.ui_button_advance_model,   'Enable', 'Off');
        set(handles.ui_button_start_auto_run,  'Enable', 'Off');
        set(handles.ui_button_reset,           'Enable', 'Off');
        set(handles.ui_edit_model_bias,        'Enable', 'Off');
        set(handles.ui_edit_nonlin_a,          'Enable', 'Off');
        set(handles.ui_edit_ens_size,          'Enable', 'Off');
        set(handles.ui_edit_inflation,         'Enable', 'Off');
        set(handles.ui_radio_button_eakf,      'Enable', 'Off');
        set(handles.ui_radio_button_enkf,      'Enable', 'Off');
        set(handles.ui_radio_button_rhf,       'Enable', 'Off');
        
    end

%% -----------------------------------------------------------------------------

    function turn_on_controls ()
        
        % Turn on all the other controls to avoid a mess
        set(handles.ui_button_advance_model,    'Enable', 'On');
        set(handles.ui_button_start_auto_run,   'Enable', 'On');
        set(handles.ui_button_reset,            'Enable', 'On');
        set(handles.ui_edit_model_bias,         'Enable', 'On');
        set(handles.ui_edit_nonlin_a,           'Enable', 'On');
        set(handles.ui_edit_ens_size,           'Enable', 'On');
        set(handles.ui_edit_inflation,          'Enable', 'On');
        set(handles.ui_radio_button_eakf,       'Enable', 'On');
        set(handles.ui_radio_button_enkf,       'Enable', 'On');
        set(handles.ui_radio_button_rhf,        'Enable', 'On');
        
    end

%% -----------------------------------------------------------------------------

    function y = calculate_rmse(x, truth)
        
        squared_error = (x - truth).^2;
        y = sqrt(mean(squared_error));
        
    end


%% -----------------------------------------------------------------------------

    function set_main_axes
        
        if (isfield(handles,'axes'))
            % 'cla reset' resets all properties of the axes except for the
            % Position and Units properties.
            cla( handles.axes,'reset');
            axes(handles.axes);
            hold off
        else
            handles.axes = axes('Units', 'Normalized', ...
                'Position', [0.050 0.382 0.333 0.400], ...
                'Color', 'White');
        end
        
        % plot some bogus items to create handles for legend
        hg_prior = plot([0 1],[0 0.1]);
        set(hg_prior, 'LineWidth', 2, 'Color', atts.green, 'Visible', 'off');
        
        hg_post  = line([0 1], [0 0.1]);
        set(hg_post, 'LineWidth', 2, 'Color', 'b', 'Visible', 'off');
        
        hg_like  = line([0 1], [0 0.1]);
        set(hg_like, 'LineWidth', 2, 'Color', 'r', 'LineStyle','--', 'Visible', 'off');
        
        legend([hg_prior hg_post hg_like],'Prior','Posterior','Likelihood');
        set(legend,'FontName', atts.fontname, 'FontSize', atts.fontsize);
        legend boxon
        
        hold on;
        axis([-10 10 -Inf Inf])

    end

%% -----------------------------------------------------------------------------

    function set_state_evolution
        
        %  axes for ensemble time series
        %  plot some items invisible just to be able to create a legend with all the
        %  potential elements.
        
        if (isfield(handles,'h_state_evolution'))
            cla( handles.h_state_evolution,'reset');
            axes(handles.h_state_evolution);
            hold off
        else
            handles.h_state_evolution = axes('Units', 'Normalized', ...
                'Position',[0.430 0.748 0.333 0.164], ...
                'Color', 'White');
        end
        
        x(1:handles.ens_size) = handles.time_step + 0.1;
        
        plot(x, handles.ens, 'b*', 'MarkerSize', 6);
        hold on;
        str = 'x_t_+_1 = x_t + (x_t+model bias) + a{\cdot}x_t{\cdot}{\mid}x_t{\mid}';
        title({str,'observation is a draw from N(0,1)'},'FontSize',18)
        
        % Include the 0 line as the truth for all times
        plot([1 100000], [0 0], 'k--');
        
        % plot the invisible stuff and capture a nice handle array for later.
        h_truth     = plot(1, 0, 'k--', 'Visible', 'off');
        h_obs       = plot(1, 0, 'r*' , 'Visible', 'off', 'MarkerSize', 10);
        h_prior     = plot(1, 0, 'g*-', 'Visible', 'off', 'MarkerSize', 6, 'Color', atts.green);
        h_posterior = plot(1, 0, 'b*' , 'Visible', 'off', 'MarkerSize', 6);
        h_evolution_handles = [h_truth h_obs h_prior h_posterior];
        
        % Want the y axis limits to take care of themselves
        set(gca, 'YLimMode', 'Auto','XTickLabel',[],'XGrid','on');
        ylabel('State','FontName', atts.fontname,'FontSize', atts.fontsize);
        
        legend(h_evolution_handles, 'Truth', 'Observation', 'Prior', 'Posterior');
        set(legend,'FontName', atts.fontname, 'FontSize', 12, ...
                   'Position',[0.821 0.770 0.118 0.148])
        legend boxon
        axis([1 10 -Inf Inf]);
        
    end

%% -----------------------------------------------------------------------------

    function set_error_spread_evolution
        % axes  for mean and spread
        % calculate rmse and spread ... expectation over a long time is
        % that they would be the same.
        
        if (isfield(handles,'h_err_spread_evolution'))
            cla( handles.h_err_spread_evolution,'reset');
            axes(handles.h_err_spread_evolution);
            hold off
        else
            handles.h_err_spread_evolution = axes('Units', 'Normalized', ...
                 'Position',[0.430 0.557 0.333 0.164], ...
                 'Color', 'White');
        end
        
        ylabel('Error, Spread','FontName', atts.fontname,'FontSize', atts.fontsize);
        axis([1 10 0 Inf]);
        set(gca,'XTickLabel',[],'XGrid','on')
        hold on;
        
    end

%% -----------------------------------------------------------------------------

    function set_kurtosis_evolution
        
        %  axes  for kurtosis
        if (isfield(handles,'h_kurtosis_evolution'))
            cla( handles.h_kurtosis_evolution,'reset');
            axes(handles.h_kurtosis_evolution);
            hold off
        else
            handles.h_kurtosis_evolution = axes('Units', 'Normalized', ...
                'Position', [0.430 0.372 0.333 0.164], ...
                'Color', 'White', ...
                'XAxisLocation','bottom');
        end
        
        axis([1 10 0 Inf]);
        ylabel('Kurtosis', 'FontName', atts.fontname, 'FontSize', atts.fontsize);
        xlabel('Timestep', 'FontName', atts.fontname, 'FontSize', atts.fontsize);
        set(gca,'XGrid', 'on')
        hold on;
        
    end

%% -------------------------------------------------------------------------

    function set_prior_histogram()
        
        %  axes for prior rank histogram
        if (isfield(handles,'h_prior_rank_histogram'))
            cla( handles.h_prior_rank_histogram,'reset');
            axes(handles.h_prior_rank_histogram);
            hold off
        else
            handles.h_prior_rank_histogram = axes('Position',[0.050 0.075 0.333 0.208]);
        end
        
        ylabel('Frequency'           ,'FontName', atts.fontname,'FontSize', atts.fontsize);
        xlabel('Rank'                ,'FontName', atts.fontname,'FontSize', atts.fontsize);
        title ('Prior Rank Histogram','FontName', atts.fontname,'FontSize', atts.fontsize);
        axis([0 handles.ens_size+2 -Inf Inf])
        set(handles.h_prior_rank_histogram,'XTick',1:(handles.ens_size+1));
        hold on;
        
    end

%% -----------------------------------------------------------------------------

    function set_posterior_histogram()
        
        %  axes for posterior rank histogram
        if (isfield(handles,'h_post_rank_histogram'))
            cla(handles.h_post_rank_histogram,'reset');
            axes(handles.h_post_rank_histogram);
            hold off
        else
            handles.h_post_rank_histogram = axes('Position',[0.43 0.075 0.333 0.208]);
        end
        
        ylabel('Frequency'               ,'FontName', atts.fontname,'FontSize', atts.fontsize);
        xlabel('Rank'                    ,'FontName', atts.fontname,'FontSize', atts.fontsize);
        title ('Posterior Rank Histogram','FontName', atts.fontname,'FontSize', atts.fontsize);
        axis([0 handles.ens_size+2 -Inf Inf])
        set(handles.h_post_rank_histogram,'XTick',1:(handles.ens_size+1));
        hold on;
        
    end

%% -----------------------------------------------------------------------------

end

% <next few lines under version control, do not edit>
% $URL$
% $Revision$
% $Date$

