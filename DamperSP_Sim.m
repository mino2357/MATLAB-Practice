classdef DamperSP_Sim  < handle
    % DamperSP_Sim 係数が線形なバネと線形なダンパーのシミュレーション
    
    properties (Constant = false, Access = public)
        % Time
        t {mustBeNumeric}
        % Spring coefficient
        k {mustBeNumeric}
        % Damper coefficient
        c {mustBeNumeric}
        % Force amplitude
        A {mustBeNumeric}
        % Wave number
        omega {mustBeNumeric}
        % const value setting
        dt       {mustBeNumeric}
        T_max    {mustBeNumeric}
        step_num {mustBeNumeric}
        % result
        time_pos {mustBeNumeric}
    end
    
    methods
        
        function ret = Spring(k, x)
            % SPRING
            ret = k * x;
        end
        
        function ret = Damper(c, vel)
            % Damper
            ret = c * vel;
        end
        
        function next_x = RK4(func, t, x, dt)
            % Classical Runge Kutta Method, Order 4.
            k1 = func(t, x);
            k2 = func(t + dt / 2.0, x + dt / 2.0 * k1);
            k3 = func(t + dt / 2.0, x + dt / 2.0 * k2);
            k4 = func(t + dt, x + dt * k3);
            next_x = x + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        end
        
        function dxdt = ode_system(t_, x_)
            % DAMPER_SP_EQ
            dxdt = zeros(2,1);
            dxdt(1) = x(2);
            dxdt(2) = - Damper(t_.c, x_(2)) - Spring(obj.k, x_(1)) + Force(obj.A, obj.omega, t_);
        end
        
        function ret = Force(A, omega, t)
            % FORCE
            ret = A * sin(omega * t);
        end
        
        function obj = DamperSP_Sim(k_, c_, A_, omega_, t0_, dt_, T_max_)
            % DAMPERSPRINGSIMULATION 物理定数の設定
            obj.k     = k_;
            obj.c     = c_;
            obj.A     = A_;
            obj.omega = omega_;
            % 時間関係の設定
            obj.dt       = dt_;
            obj.T_max    = T_max_;
            obj.step_num = T_max_ / dt_;
            obj.t        = t0_;
            obj.time_pos = zeros(2, T_max_ / dt_);
        end
        
        function calc_ode(obj, x0, v0)
            % calc_ode 
            y1 = [x0; v0];
            
            for i = 1:obj.step_num
                obj.time_pos(1,i) = obj.t;
                obj.time_pos(2,i) = y1(1);
                y2 = RK4(@ode_system, obj.t, y1, obj.dt);
                obj.t  = obj.t + obj.dt;
                y1 = y2;
            end
        end
        
        function plot_result(obj)
            % plot
            plot(obj.time_pos(1,:), obj.time_pos(2,:));
        end
    end
end

