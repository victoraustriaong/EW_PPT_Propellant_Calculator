classdef uPPT_Performance_Simulator_V5_exported_uPPT_params < matlab.apps.AppBase
%This is GUI developed to simulate performance of uPPT. Upon running this
%matlab code, the values on the fields are pre-populated by the parameters
%of our uPPT. Moreover, input fields can be modified to conduct more
%simulations.
%Last run: 4/4/2020
%Author: Victor Ong

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        ThrusterSpecificationsPanel     matlab.ui.container.Panel
        GapHeightmLabel                 matlab.ui.control.Label
        GapHeightmEditField             matlab.ui.control.NumericEditField
        PlateWidthmEditFieldLabel       matlab.ui.control.Label
        PlateWidthmEditField            matlab.ui.control.NumericEditField
        ParallelPlateGeometryLabel      matlab.ui.control.Label
        PlateThickmEditFieldLabel       matlab.ui.control.Label
        PlateThickmEditField            matlab.ui.control.NumericEditField
        ChannelLengthmEditFieldLabel    matlab.ui.control.Label
        ChannelLengthmEditField         matlab.ui.control.NumericEditField
        CircuitParametersLabel          matlab.ui.control.Label
        CapacitanceFEditFieldLabel      matlab.ui.control.Label
        CapacitanceFEditField           matlab.ui.control.NumericEditField
        InitialVoltageVEditFieldLabel   matlab.ui.control.Label
        InitialVoltageVEditField        matlab.ui.control.NumericEditField
        CapacitorLcHEditFieldLabel      matlab.ui.control.Label
        CapacitorLcHEditField           matlab.ui.control.NumericEditField
        InductanceHenryLabel            matlab.ui.control.Label
        LeadWireLeHEditFieldLabel       matlab.ui.control.Label
        LeadWireLeHEditField            matlab.ui.control.NumericEditField
        ElectrodeLpeHEditFieldLabel     matlab.ui.control.Label
        ElectrodeLpeHEditField          matlab.ui.control.NumericEditField
        ResistanceOhmsLabel             matlab.ui.control.Label
        CapacitorRcEditFieldLabel       matlab.ui.control.Label
        CapacitorRcEditField            matlab.ui.control.NumericEditField
        LeadWireReEditFieldLabel        matlab.ui.control.Label
        LeadWireReEditField             matlab.ui.control.NumericEditField
        ElectrodeRpeEditFieldLabel      matlab.ui.control.Label
        ElectrodeRpeEditField           matlab.ui.control.NumericEditField
        PlasmaRpEditFieldLabel          matlab.ui.control.Label
        PlasmaRpEditField               matlab.ui.control.NumericEditField
        PlasmaResistanceModelLabel      matlab.ui.control.Label
        UsePlasmaResistanceModelCheckBox  matlab.ui.control.CheckBox
        CharacteristicPulseTimesEditFieldLabel  matlab.ui.control.Label
        CharacteristicPulseTimesEditField  matlab.ui.control.NumericEditField
        ElectronNumDensity1m3Label      matlab.ui.control.Label
        ElectronNumDensity1m3EditField  matlab.ui.control.NumericEditField
        PlasmaElectronTempeVEditFieldLabel  matlab.ui.control.Label
        PlasmaElectronTempeVEditField   matlab.ui.control.NumericEditField
        MicroPulsedPlasmaThrusterSimulationLabel  matlab.ui.control.Label
        ModelSpecificationsPanel        matlab.ui.container.Panel
        InitialMasst0kgEditFieldLabel   matlab.ui.control.Label
        InitialMasst0kgEditField        matlab.ui.control.NumericEditField
        EntrainedMasskgEditFieldLabel   matlab.ui.control.Label
        EntrainedMasskgEditField        matlab.ui.control.NumericEditField
        MassdistloadingparamEditFieldLabel  matlab.ui.control.Label
        MassdistloadingparamEditField   matlab.ui.control.NumericEditField
        ConstantMassSlugOperationCheckBox  matlab.ui.control.CheckBox
        SolverSpecificationsPanel       matlab.ui.container.Panel
        FinalTimesEditFieldLabel        matlab.ui.control.Label
        FinalTimesEditField             matlab.ui.control.NumericEditField
        TimestepsEditFieldLabel         matlab.ui.control.Label
        TimestepsEditField              matlab.ui.control.NumericEditField
        ToleranceEditFieldLabel         matlab.ui.control.Label
        ToleranceEditField              matlab.ui.control.NumericEditField
        SolverDropDownLabel             matlab.ui.control.Label
        SolverDropDown                  matlab.ui.control.DropDown
        SimulationResultsPanel          matlab.ui.container.Panel
        ExhaustVelocitymsEditFieldLabel  matlab.ui.control.Label
        ExhaustVelocitymsEditField      matlab.ui.control.NumericEditField
        SpecificImpulsesEditFieldLabel  matlab.ui.control.Label
        SpecificImpulsesEditField       matlab.ui.control.NumericEditField
        ImpulseBituNsEditFieldLabel     matlab.ui.control.Label
        ImpulseBituNsEditField          matlab.ui.control.NumericEditField
        EfficiencyEditFieldLabel        matlab.ui.control.Label
        EfficiencyEditField             matlab.ui.control.NumericEditField
        ThrustPowerEditFieldLabel       matlab.ui.control.Label
        ThrustPowerEditField            matlab.ui.control.NumericEditField
        RunSimulationButton             matlab.ui.control.Button
        CommentsLabel                   matlab.ui.control.Label
        CommentBox                      matlab.ui.control.TextArea
        PlotsPanel                      matlab.ui.container.Panel
        UIAxes_VoltCur                  matlab.ui.control.UIAxes
        UIAxes_Energy                   matlab.ui.control.UIAxes
        UIAxes_PosVel                   matlab.ui.control.UIAxes
    end

    
    properties (Access = public)
        %Allows variable to be seen by other GUI functions
        V0;C;Rc;Re;Rpe;Rp;R_tot;Lc;Le;Lpe;Te;ne;tau;m0;m_ent;m_t;alpha;ch_length;h;w;t_el;constantMassModelCheck;plasmaModelCheck;
        %Constant        
        u0 = 1.2566e-6; %Wb/Am
    
    end
    
    methods (Access=public)
        %Function needed by the ODE solver
        function xdot=func_xdot(app,t,x)
                       
            if app.constantMassModelCheck==1
                app.m_t = app.m0;
            else
                app.m_t = app.m0 + app.m_ent * (1 - (1- x(1)/app.ch_length)^(1/(1-app.alpha)));
            end
            
            %Inductance which includes electrode thickness
           app.Lpe = (app.u0/pi) * (1.5 + log (app.h/(app.w+app.t_el)));
            
            %Inductance which does not includes electrode thickness
            %app.Lpe = app.u0 *(app.h/app.w);
            Lpe_t = app.Lpe*x(1);
            
            L_tot_t = app.Lc + app.Le + Lpe_t;
            
                        
            if app.plasmaModelCheck==1
                app.Rp = ((8.08 * app.h) /(app.w*(app.Te^0.75)))  * sqrt((app.u0 * log (1.24e7 * ((app.Te^3) / app.ne )^0.5 ))/app.tau);
            end
            
            app.R_tot = app.Rc + app.Re + app.Rpe + app.Rp;
            
            %Solution:
            xdot(1) = x(3); %Current Sheet Position
            xdot(2) = x(4); %Charge on capacitor
            xdot(3) = (0.5/app.m_t)*(app.u0*app.h/app.w)*(x(4))^2; %Velocity     
            xdot(4) = (app.V0 - (x(2)/app.C) - (x(4)*app.R_tot) - (app.Lpe*x(3)*x(4)))/L_tot_t; %Current
            
            xdot(5) = app.R_tot*x(4)^2; %Ohmic energy
            xdot(6) = 0.5*app.u0*(app.h/app.w)*x(4)^2;%Impulse bit
            
            xdot = xdot';
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: SolverDropDown
        function SolverDropDownValueChanged(app, event)
           
        end

        % Button pushed function: RunSimulationButton
        function RunSimulationButtonPushed(app, event)
            
            
            
            solverVal = app.SolverDropDown.Value;
            finalTime = app.FinalTimesEditField.Value;
            timesteps = app.TimestepsEditField.Value;
            tolerance = app.ToleranceEditField.Value;

            
            app.V0 = app.InitialVoltageVEditField.Value;
            app.C = app.CapacitanceFEditField.Value;
            
            %Resistance:
            app.Rc = app.CapacitorRcEditField.Value;
            app.Re = app.LeadWireReEditField.Value;
            app.Rpe = app.ElectrodeRpeEditField.Value;
            app.Rp = app.PlasmaRpEditField.Value;
            
            %Inductance:
            app.Lc = app.CapacitorLcHEditField.Value;
            app.Le = app.LeadWireLeHEditField.Value;
            app.Lpe = app.ElectrodeLpeHEditField.Value;
            
            %Plasma model parameters:
            app.Te = (app.PlasmaElectronTempeVEditField.Value * 1.60218e-19)/1.3807e-23 ; %Kelvin;
            app.ne = app.ElectronNumDensity1m3EditField.Value;
            app.tau = app.CharacteristicPulseTimesEditField.Value;
            
            %Distributed mass model parameters:
            app.m0 = app.InitialMasst0kgEditField.Value;
            app.m_ent = app.EntrainedMasskgEditField.Value;
            app.alpha = app.MassdistloadingparamEditField.Value;
            app.ch_length = app.ChannelLengthmEditField.Value;
                                
            %Geometry:
            app.h = app.GapHeightmEditField.Value;
            app.w = app.PlateWidthmEditField.Value;
            app.t_el = app.PlateThickmEditField.Value;
            
            
            app.constantMassModelCheck = app.ConstantMassSlugOperationCheckBox.Value;
            
            app.plasmaModelCheck = app.UsePlasmaResistanceModelCheckBox.Value;
                        
            
            options = odeset('RelTol',tolerance,'AbsTol', tolerance);           
            tt=linspace(0,finalTime,timesteps);         
                     
            if solverVal=="ODE45"
                [t,x] = ode45(@(t,x) func_xdot(app,t,x),tt, [0 0 0 0 0 0],options); 
            elseif solverVal=="ODE23S"
                [t,x] = ode23s(@(t,x) func_xdot(app,t,x),tt, [0 0 0 0],options); 
            elseif solverVal=="ODE15S"
                [t,x] = ode15s(@(t,x) func_xdot(app,t,x),tt, [0 0 0 0],options); 
            end
            
            %Results as function of time
            position = x(:,1);
            voltage = app.V0-(x(:,2)./app.C);
            velocity = x(:,3);
            current = x(:,4);
            
            
            %Performance calculations:
            channelExitTime = interp1(position,t,app.ch_length);
            channelExitVelocity = interp1(t,velocity,channelExitTime);
            iBit = interp1(t,x(:,6),channelExitTime)/(1e-6);                       
            
            if isnan(channelExitVelocity)
                app.CommentBox.Value = "The result value with respect to the exit time at channel length is NaN. Try increasing the simulation time.";
                
                app.ExhaustVelocitymsEditField.Value = 0;
                app.SpecificImpulsesEditField.Value = 0;
                app.EfficiencyEditField.Value = 0;
                app.ImpulseBituNsEditField.Value = 0;
            else
                app.ExhaustVelocitymsEditField.Value = channelExitVelocity;%exVel;
                app.SpecificImpulsesEditField.Value = channelExitVelocity/9.81;
                app.EfficiencyEditField.Value = (app.m_t*(channelExitVelocity^2) / (app.C*app.V0^2))*100;
                app.ImpulseBituNsEditField.Value = iBit;%(channelExitVelocity * app.m_t)/(1e-6);  
                app.CommentBox.Value = "";
            end
            
            
            %Plots
            
            yyaxis(app.UIAxes_PosVel,'left')
            plot(app.UIAxes_PosVel,t,position,channelExitTime, app.ch_length,'o')%Position
            app.UIAxes_PosVel.YLabel.String='Position (m)';
            yyaxis(app.UIAxes_PosVel,'right')
            plot(app.UIAxes_PosVel,t,velocity, channelExitTime,channelExitVelocity,'o')%Velocity
            app.UIAxes_PosVel.YLabel.String='Velocity (m/s)';
            legend(app.UIAxes_PosVel,[{'Position','Exit Time'},{'Velocity.','Exit Time'}])
            
            
            yyaxis(app.UIAxes_VoltCur,'left')
            plot(app.UIAxes_VoltCur,t,voltage,'g')%Voltage
            app.UIAxes_VoltCur.YLabel.String='Voltage (V)';
            yyaxis(app.UIAxes_VoltCur,'right')
            plot(app.UIAxes_VoltCur,t,current)%Current
            app.UIAxes_VoltCur.YLabel.String='Current (A)';
            legend(app.UIAxes_VoltCur,[{'Voltage'},{'Current'}])
            
            %Energy
            cap_energy = (0.5*app.C)* (voltage.^2); %Capacitor energy
            mag_energy = 0.5 * (app.Lc + app.Le + (app.u0*app.h/app.w) .* position) .* (current.^2);
            kin_energy = 0.5* app.m0 *(velocity.^2);
            ohmic_energy = x(:,5);
            total_energy = cap_energy + mag_energy + kin_energy + ohmic_energy;
                                  
            intcap = interp1(t,cap_energy,channelExitTime);
            intcap0 = interp1(t,cap_energy,0);
            
            TP = iBit / intcap0;
            app.ThrustPowerEditField.Value = TP;
            
%             intmag = interp1(t,mag_energy,channelExitTime);
%             intkin = interp1(t,kin_energy,channelExitTime);                 
                         
            plot(app.UIAxes_Energy,t,cap_energy,t,mag_energy,'r',t,kin_energy, 'g', t,ohmic_energy, t, total_energy)
            app.UIAxes_Energy.YLabel.String='Energy (J)';
            legend(app.UIAxes_Energy,'Capacitor Energy', 'Magnetic Field Induced Energy','Kinetic Energy', 'Ohmic Energy', 'Total Energy')
                        
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1168 774];
            app.UIFigure.Name = 'UI Figure';

            % Create ThrusterSpecificationsPanel
            app.ThrusterSpecificationsPanel = uipanel(app.UIFigure);
            app.ThrusterSpecificationsPanel.Title = '1.) Thruster Specifications';
            app.ThrusterSpecificationsPanel.FontWeight = 'bold';
            app.ThrusterSpecificationsPanel.FontSize = 14;
            app.ThrusterSpecificationsPanel.Position = [9 14 253 715];

            % Create GapHeightmLabel
            app.GapHeightmLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.GapHeightmLabel.HorizontalAlignment = 'right';
            app.GapHeightmLabel.Position = [3 639 91 22];
            app.GapHeightmLabel.Text = 'Gap Height (m):';

            % Create GapHeightmEditField
            app.GapHeightmEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.GapHeightmEditField.Position = [130 639 68 22];
            app.GapHeightmEditField.Value = 0.02;

            % Create PlateWidthmEditFieldLabel
            app.PlateWidthmEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.PlateWidthmEditFieldLabel.HorizontalAlignment = 'right';
            app.PlateWidthmEditFieldLabel.Position = [3 609 92 22];
            app.PlateWidthmEditFieldLabel.Text = 'Plate Width (m):';

            % Create PlateWidthmEditField
            app.PlateWidthmEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.PlateWidthmEditField.Position = [130 609 68 22];
            app.PlateWidthmEditField.Value = 0.01;

            % Create ParallelPlateGeometryLabel
            app.ParallelPlateGeometryLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.ParallelPlateGeometryLabel.FontSize = 13;
            app.ParallelPlateGeometryLabel.FontWeight = 'bold';
            app.ParallelPlateGeometryLabel.Position = [5 663 155 22];
            app.ParallelPlateGeometryLabel.Text = 'Parallel Plate Geometry:';

            % Create PlateThickmEditFieldLabel
            app.PlateThickmEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.PlateThickmEditFieldLabel.HorizontalAlignment = 'right';
            app.PlateThickmEditFieldLabel.Position = [3 578 93 22];
            app.PlateThickmEditFieldLabel.Text = 'Plate Thick. (m):';

            % Create PlateThickmEditField
            app.PlateThickmEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.PlateThickmEditField.Position = [130 577 68 22];
            app.PlateThickmEditField.Value = 0.002;

            % Create ChannelLengthmEditFieldLabel
            app.ChannelLengthmEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.ChannelLengthmEditFieldLabel.HorizontalAlignment = 'right';
            app.ChannelLengthmEditFieldLabel.Position = [1 546 115 22];
            app.ChannelLengthmEditFieldLabel.Text = 'Channel Length (m):';

            % Create ChannelLengthmEditField
            app.ChannelLengthmEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.ChannelLengthmEditField.Position = [130 546 68 22];
            app.ChannelLengthmEditField.Value = 0.021;

            % Create CircuitParametersLabel
            app.CircuitParametersLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.CircuitParametersLabel.FontSize = 13;
            app.CircuitParametersLabel.FontWeight = 'bold';
            app.CircuitParametersLabel.Position = [5 506 125 22];
            app.CircuitParametersLabel.Text = 'Circuit Parameters:';

            % Create CapacitanceFEditFieldLabel
            app.CapacitanceFEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.CapacitanceFEditFieldLabel.HorizontalAlignment = 'right';
            app.CapacitanceFEditFieldLabel.Position = [24 452 94 22];
            app.CapacitanceFEditFieldLabel.Text = 'Capacitance (F):';

            % Create CapacitanceFEditField
            app.CapacitanceFEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.CapacitanceFEditField.Position = [159 452 67 22];
            app.CapacitanceFEditField.Value = 3e-06;

            % Create InitialVoltageVEditFieldLabel
            app.InitialVoltageVEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.InitialVoltageVEditFieldLabel.HorizontalAlignment = 'right';
            app.InitialVoltageVEditFieldLabel.Position = [24 482 100 22];
            app.InitialVoltageVEditFieldLabel.Text = 'Initial Voltage (V):';

            % Create InitialVoltageVEditField
            app.InitialVoltageVEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.InitialVoltageVEditField.Position = [159 482 67 22];
            app.InitialVoltageVEditField.Value = 1500;

            % Create CapacitorLcHEditFieldLabel
            app.CapacitorLcHEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.CapacitorLcHEditFieldLabel.HorizontalAlignment = 'right';
            app.CapacitorLcHEditFieldLabel.Position = [24 402 96 22];
            app.CapacitorLcHEditFieldLabel.Text = 'Capacitor Lc (H):';

            % Create CapacitorLcHEditField
            app.CapacitorLcHEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.CapacitorLcHEditField.Position = [160 402 67 22];
            app.CapacitorLcHEditField.Value = 3.5e-08;

            % Create InductanceHenryLabel
            app.InductanceHenryLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.InductanceHenryLabel.FontWeight = 'bold';
            app.InductanceHenryLabel.Position = [17 423 118 22];
            app.InductanceHenryLabel.Text = 'Inductance (Henry):';

            % Create LeadWireLeHEditFieldLabel
            app.LeadWireLeHEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.LeadWireLeHEditFieldLabel.HorizontalAlignment = 'right';
            app.LeadWireLeHEditFieldLabel.Position = [24 371 112 22];
            app.LeadWireLeHEditFieldLabel.Text = 'Lead & Wire Le (H):';

            % Create LeadWireLeHEditField
            app.LeadWireLeHEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.LeadWireLeHEditField.Position = [159 371 67 22];

            % Create ElectrodeLpeHEditFieldLabel
            app.ElectrodeLpeHEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.ElectrodeLpeHEditFieldLabel.HorizontalAlignment = 'right';
            app.ElectrodeLpeHEditFieldLabel.Enable = 'off';
            app.ElectrodeLpeHEditFieldLabel.Visible = 'off';
            app.ElectrodeLpeHEditFieldLabel.Position = [24 340 103 22];
            app.ElectrodeLpeHEditFieldLabel.Text = 'Electrode Lpe (H):';

            % Create ElectrodeLpeHEditField
            app.ElectrodeLpeHEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.ElectrodeLpeHEditField.Enable = 'off';
            app.ElectrodeLpeHEditField.Visible = 'off';
            app.ElectrodeLpeHEditField.Position = [159 340 67 22];

            % Create ResistanceOhmsLabel
            app.ResistanceOhmsLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.ResistanceOhmsLabel.FontWeight = 'bold';
            app.ResistanceOhmsLabel.Position = [17 309 118 22];
            app.ResistanceOhmsLabel.Text = 'Resistance (Ohms):';

            % Create CapacitorRcEditFieldLabel
            app.CapacitorRcEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.CapacitorRcEditFieldLabel.HorizontalAlignment = 'right';
            app.CapacitorRcEditFieldLabel.Position = [24 286 78 22];
            app.CapacitorRcEditFieldLabel.Text = 'Capacitor Rc:';

            % Create CapacitorRcEditField
            app.CapacitorRcEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.CapacitorRcEditField.Position = [159 286 67 22];
            app.CapacitorRcEditField.Value = 0.03;

            % Create LeadWireReEditFieldLabel
            app.LeadWireReEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.LeadWireReEditFieldLabel.HorizontalAlignment = 'right';
            app.LeadWireReEditFieldLabel.Position = [24 255 94 22];
            app.LeadWireReEditFieldLabel.Text = 'Lead & Wire Re:';

            % Create LeadWireReEditField
            app.LeadWireReEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.LeadWireReEditField.Position = [159 255 67 22];

            % Create ElectrodeRpeEditFieldLabel
            app.ElectrodeRpeEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.ElectrodeRpeEditFieldLabel.HorizontalAlignment = 'right';
            app.ElectrodeRpeEditFieldLabel.Position = [24 224 85 22];
            app.ElectrodeRpeEditFieldLabel.Text = 'Electrode Rpe:';

            % Create ElectrodeRpeEditField
            app.ElectrodeRpeEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.ElectrodeRpeEditField.Position = [159 224 67 22];

            % Create PlasmaRpEditFieldLabel
            app.PlasmaRpEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.PlasmaRpEditFieldLabel.HorizontalAlignment = 'right';
            app.PlasmaRpEditFieldLabel.Position = [24 193 68 22];
            app.PlasmaRpEditFieldLabel.Text = 'Plasma Rp:';

            % Create PlasmaRpEditField
            app.PlasmaRpEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.PlasmaRpEditField.ValueDisplayFormat = '%.7f';
            app.PlasmaRpEditField.Position = [160 193 67 22];
            app.PlasmaRpEditField.Value = 0.0130392;

            % Create PlasmaResistanceModelLabel
            app.PlasmaResistanceModelLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.PlasmaResistanceModelLabel.FontSize = 13;
            app.PlasmaResistanceModelLabel.FontWeight = 'bold';
            app.PlasmaResistanceModelLabel.Position = [3 159 169 22];
            app.PlasmaResistanceModelLabel.Text = 'Plasma Resistance Model:';

            % Create UsePlasmaResistanceModelCheckBox
            app.UsePlasmaResistanceModelCheckBox = uicheckbox(app.ThrusterSpecificationsPanel);
            app.UsePlasmaResistanceModelCheckBox.Text = 'Use Plasma Resistance Model';
            app.UsePlasmaResistanceModelCheckBox.Position = [9 138 185 22];
            app.UsePlasmaResistanceModelCheckBox.Value = true;

            % Create CharacteristicPulseTimesEditFieldLabel
            app.CharacteristicPulseTimesEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.CharacteristicPulseTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.CharacteristicPulseTimesEditFieldLabel.Position = [3 98 162 22];
            app.CharacteristicPulseTimesEditFieldLabel.Text = 'Characteristic Pulse Time (s):';

            % Create CharacteristicPulseTimesEditField
            app.CharacteristicPulseTimesEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.CharacteristicPulseTimesEditField.Position = [180 98 67 22];
            app.CharacteristicPulseTimesEditField.Value = 4e-07;

            % Create ElectronNumDensity1m3Label
            app.ElectronNumDensity1m3Label = uilabel(app.ThrusterSpecificationsPanel);
            app.ElectronNumDensity1m3Label.HorizontalAlignment = 'right';
            app.ElectronNumDensity1m3Label.Position = [3 69 169 22];
            app.ElectronNumDensity1m3Label.Text = 'Electron Num Density (1/m^3):';

            % Create ElectronNumDensity1m3EditField
            app.ElectronNumDensity1m3EditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.ElectronNumDensity1m3EditField.Position = [180 69 67 22];
            app.ElectronNumDensity1m3EditField.Value = 1e+21;

            % Create PlasmaElectronTempeVEditFieldLabel
            app.PlasmaElectronTempeVEditFieldLabel = uilabel(app.ThrusterSpecificationsPanel);
            app.PlasmaElectronTempeVEditFieldLabel.HorizontalAlignment = 'right';
            app.PlasmaElectronTempeVEditFieldLabel.Position = [4 37 158 22];
            app.PlasmaElectronTempeVEditFieldLabel.Text = 'Plasma Electron Temp. (eV):';

            % Create PlasmaElectronTempeVEditField
            app.PlasmaElectronTempeVEditField = uieditfield(app.ThrusterSpecificationsPanel, 'numeric');
            app.PlasmaElectronTempeVEditField.Position = [180 37 67 22];
            app.PlasmaElectronTempeVEditField.Value = 1.5;

            % Create MicroPulsedPlasmaThrusterSimulationLabel
            app.MicroPulsedPlasmaThrusterSimulationLabel = uilabel(app.UIFigure);
            app.MicroPulsedPlasmaThrusterSimulationLabel.FontSize = 16;
            app.MicroPulsedPlasmaThrusterSimulationLabel.FontWeight = 'bold';
            app.MicroPulsedPlasmaThrusterSimulationLabel.Position = [9 743 325 22];
            app.MicroPulsedPlasmaThrusterSimulationLabel.Text = 'Micro Pulsed Plasma Thruster Simulation';

            % Create ModelSpecificationsPanel
            app.ModelSpecificationsPanel = uipanel(app.UIFigure);
            app.ModelSpecificationsPanel.Title = '2.) Model Specifications';
            app.ModelSpecificationsPanel.FontWeight = 'bold';
            app.ModelSpecificationsPanel.FontSize = 14;
            app.ModelSpecificationsPanel.Position = [269 571 224 158];

            % Create InitialMasst0kgEditFieldLabel
            app.InitialMasst0kgEditFieldLabel = uilabel(app.ModelSpecificationsPanel);
            app.InitialMasst0kgEditFieldLabel.HorizontalAlignment = 'right';
            app.InitialMasst0kgEditFieldLabel.Position = [1 80 117 22];
            app.InitialMasst0kgEditFieldLabel.Text = 'Initial Mass, t=0 (kg):';

            % Create InitialMasst0kgEditField
            app.InitialMasst0kgEditField = uieditfield(app.ModelSpecificationsPanel, 'numeric');
            app.InitialMasst0kgEditField.Position = [150 80 68 22];
            app.InitialMasst0kgEditField.Value = 1e-08;

            % Create EntrainedMasskgEditFieldLabel
            app.EntrainedMasskgEditFieldLabel = uilabel(app.ModelSpecificationsPanel);
            app.EntrainedMasskgEditFieldLabel.HorizontalAlignment = 'right';
            app.EntrainedMasskgEditFieldLabel.Position = [2 51 116 22];
            app.EntrainedMasskgEditFieldLabel.Text = 'Entrained Mass (kg):';

            % Create EntrainedMasskgEditField
            app.EntrainedMasskgEditField = uieditfield(app.ModelSpecificationsPanel, 'numeric');
            app.EntrainedMasskgEditField.Position = [150 51 68 22];

            % Create MassdistloadingparamEditFieldLabel
            app.MassdistloadingparamEditFieldLabel = uilabel(app.ModelSpecificationsPanel);
            app.MassdistloadingparamEditFieldLabel.HorizontalAlignment = 'right';
            app.MassdistloadingparamEditFieldLabel.Position = [1 21 142 22];
            app.MassdistloadingparamEditFieldLabel.Text = 'Mass dist. loading param:';

            % Create MassdistloadingparamEditField
            app.MassdistloadingparamEditField = uieditfield(app.ModelSpecificationsPanel, 'numeric');
            app.MassdistloadingparamEditField.ValueDisplayFormat = '%.5f';
            app.MassdistloadingparamEditField.Position = [150 21 68 22];
            app.MassdistloadingparamEditField.Value = 0.99999;

            % Create ConstantMassSlugOperationCheckBox
            app.ConstantMassSlugOperationCheckBox = uicheckbox(app.ModelSpecificationsPanel);
            app.ConstantMassSlugOperationCheckBox.Text = 'Constant Mass (Slug Operation) ';
            app.ConstantMassSlugOperationCheckBox.Position = [7 106 198 22];
            app.ConstantMassSlugOperationCheckBox.Value = true;

            % Create SolverSpecificationsPanel
            app.SolverSpecificationsPanel = uipanel(app.UIFigure);
            app.SolverSpecificationsPanel.Title = '3.) Solver Specifications';
            app.SolverSpecificationsPanel.FontWeight = 'bold';
            app.SolverSpecificationsPanel.FontSize = 14;
            app.SolverSpecificationsPanel.Position = [270 385 224 179];

            % Create FinalTimesEditFieldLabel
            app.FinalTimesEditFieldLabel = uilabel(app.SolverSpecificationsPanel);
            app.FinalTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.FinalTimesEditFieldLabel.Position = [9 88 82 22];
            app.FinalTimesEditFieldLabel.Text = 'Final Time (s):';

            % Create FinalTimesEditField
            app.FinalTimesEditField = uieditfield(app.SolverSpecificationsPanel, 'numeric');
            app.FinalTimesEditField.Position = [148 88 68 22];
            app.FinalTimesEditField.Value = 1.8e-05;

            % Create TimestepsEditFieldLabel
            app.TimestepsEditFieldLabel = uilabel(app.SolverSpecificationsPanel);
            app.TimestepsEditFieldLabel.HorizontalAlignment = 'right';
            app.TimestepsEditFieldLabel.Position = [9 59 64 22];
            app.TimestepsEditFieldLabel.Text = 'Timesteps:';

            % Create TimestepsEditField
            app.TimestepsEditField = uieditfield(app.SolverSpecificationsPanel, 'numeric');
            app.TimestepsEditField.Position = [148 59 68 22];
            app.TimestepsEditField.Value = 1000;

            % Create ToleranceEditFieldLabel
            app.ToleranceEditFieldLabel = uilabel(app.SolverSpecificationsPanel);
            app.ToleranceEditFieldLabel.HorizontalAlignment = 'right';
            app.ToleranceEditFieldLabel.Position = [9 29 61 22];
            app.ToleranceEditFieldLabel.Text = 'Tolerance:';

            % Create ToleranceEditField
            app.ToleranceEditField = uieditfield(app.SolverSpecificationsPanel, 'numeric');
            app.ToleranceEditField.Position = [148 29 68 22];
            app.ToleranceEditField.Value = 1e-06;

            % Create SolverDropDownLabel
            app.SolverDropDownLabel = uilabel(app.SolverSpecificationsPanel);
            app.SolverDropDownLabel.HorizontalAlignment = 'right';
            app.SolverDropDownLabel.Position = [1 125 43 22];
            app.SolverDropDownLabel.Text = 'Solver:';

            % Create SolverDropDown
            app.SolverDropDown = uidropdown(app.SolverSpecificationsPanel);
            app.SolverDropDown.Items = {'ODE45', 'ODE23S', 'ODE15S'};
            app.SolverDropDown.ValueChangedFcn = createCallbackFcn(app, @SolverDropDownValueChanged, true);
            app.SolverDropDown.Position = [69 125 147 22];
            app.SolverDropDown.Value = 'ODE45';

            % Create SimulationResultsPanel
            app.SimulationResultsPanel = uipanel(app.UIFigure);
            app.SimulationResultsPanel.Title = '4.) Simulation Results';
            app.SimulationResultsPanel.FontWeight = 'bold';
            app.SimulationResultsPanel.FontSize = 14;
            app.SimulationResultsPanel.Position = [271 183 224 193];

            % Create ExhaustVelocitymsEditFieldLabel
            app.ExhaustVelocitymsEditFieldLabel = uilabel(app.SimulationResultsPanel);
            app.ExhaustVelocitymsEditFieldLabel.HorizontalAlignment = 'right';
            app.ExhaustVelocitymsEditFieldLabel.Position = [5 138 128 22];
            app.ExhaustVelocitymsEditFieldLabel.Text = 'Exhaust Velocity (m/s):';

            % Create ExhaustVelocitymsEditField
            app.ExhaustVelocitymsEditField = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.ExhaustVelocitymsEditField.ValueDisplayFormat = '%.3f';
            app.ExhaustVelocitymsEditField.Position = [149 138 68 22];

            % Create SpecificImpulsesEditFieldLabel
            app.SpecificImpulsesEditFieldLabel = uilabel(app.SimulationResultsPanel);
            app.SpecificImpulsesEditFieldLabel.HorizontalAlignment = 'right';
            app.SpecificImpulsesEditFieldLabel.Position = [5 109 114 22];
            app.SpecificImpulsesEditFieldLabel.Text = 'Specific Impulse (s):';

            % Create SpecificImpulsesEditField
            app.SpecificImpulsesEditField = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.SpecificImpulsesEditField.ValueDisplayFormat = '%.3f';
            app.SpecificImpulsesEditField.Position = [150 109 68 22];

            % Create ImpulseBituNsEditFieldLabel
            app.ImpulseBituNsEditFieldLabel = uilabel(app.SimulationResultsPanel);
            app.ImpulseBituNsEditFieldLabel.HorizontalAlignment = 'right';
            app.ImpulseBituNsEditFieldLabel.Position = [5 79 105 22];
            app.ImpulseBituNsEditFieldLabel.Text = 'Impulse Bit (uN-s):';

            % Create ImpulseBituNsEditField
            app.ImpulseBituNsEditField = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.ImpulseBituNsEditField.ValueDisplayFormat = '%.3f';
            app.ImpulseBituNsEditField.Position = [150 79 68 22];

            % Create EfficiencyEditFieldLabel
            app.EfficiencyEditFieldLabel = uilabel(app.SimulationResultsPanel);
            app.EfficiencyEditFieldLabel.HorizontalAlignment = 'right';
            app.EfficiencyEditFieldLabel.Position = [6 48 82 22];
            app.EfficiencyEditFieldLabel.Text = 'Efficiency (%):';

            % Create EfficiencyEditField
            app.EfficiencyEditField = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.EfficiencyEditField.ValueDisplayFormat = '%.2f';
            app.EfficiencyEditField.Position = [150 48 68 22];

            % Create ThrustPowerEditFieldLabel
            app.ThrustPowerEditFieldLabel = uilabel(app.SimulationResultsPanel);
            app.ThrustPowerEditFieldLabel.HorizontalAlignment = 'right';
            app.ThrustPowerEditFieldLabel.Position = [6 16 80 22];
            app.ThrustPowerEditFieldLabel.Text = 'Thrust/Power:';

            % Create ThrustPowerEditField
            app.ThrustPowerEditField = uieditfield(app.SimulationResultsPanel, 'numeric');
            app.ThrustPowerEditField.ValueDisplayFormat = '%.2f';
            app.ThrustPowerEditField.Position = [151 16 68 22];

            % Create RunSimulationButton
            app.RunSimulationButton = uibutton(app.UIFigure, 'push');
            app.RunSimulationButton.ButtonPushedFcn = createCallbackFcn(app, @RunSimulationButtonPushed, true);
            app.RunSimulationButton.FontSize = 18;
            app.RunSimulationButton.FontWeight = 'bold';
            app.RunSimulationButton.Position = [292 31 182 29];
            app.RunSimulationButton.Text = 'Run Simulation';

            % Create CommentsLabel
            app.CommentsLabel = uilabel(app.UIFigure);
            app.CommentsLabel.FontWeight = 'bold';
            app.CommentsLabel.Position = [271 150 72 22];
            app.CommentsLabel.Text = 'Comments:';

            % Create CommentBox
            app.CommentBox = uitextarea(app.UIFigure);
            app.CommentBox.Editable = 'off';
            app.CommentBox.Position = [273 89 222 60];

            % Create PlotsPanel
            app.PlotsPanel = uipanel(app.UIFigure);
            app.PlotsPanel.Title = 'Plots';
            app.PlotsPanel.FontWeight = 'bold';
            app.PlotsPanel.FontSize = 14;
            app.PlotsPanel.Position = [503 14 658 751];

            % Create UIAxes_VoltCur
            app.UIAxes_VoltCur = uiaxes(app.PlotsPanel);
            title(app.UIAxes_VoltCur, 'Voltage and Current vs. Time')
            xlabel(app.UIAxes_VoltCur, 'Time (s)')
            ylabel(app.UIAxes_VoltCur, '')
            app.UIAxes_VoltCur.Position = [10 255 629 223];

            % Create UIAxes_Energy
            app.UIAxes_Energy = uiaxes(app.PlotsPanel);
            title(app.UIAxes_Energy, 'Energy vs. Time')
            xlabel(app.UIAxes_Energy, 'Time (s)')
            ylabel(app.UIAxes_Energy, '')
            app.UIAxes_Energy.PlotBoxAspectRatio = [3.86153846153846 1 1];
            app.UIAxes_Energy.Position = [10 7 629 232];

            % Create UIAxes_PosVel
            app.UIAxes_PosVel = uiaxes(app.PlotsPanel);
            title(app.UIAxes_PosVel, 'Position and Velocity vs. Time')
            xlabel(app.UIAxes_PosVel, 'Time (s)')
            ylabel(app.UIAxes_PosVel, '')
            app.UIAxes_PosVel.Position = [10 492 629 223];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = uPPT_Performance_Simulator_V5_exported_uPPT_params

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end