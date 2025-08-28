classdef createTWPA
    properties
        fsim {mustBeNumeric}
        ksim {mustBeNumeric}
        gsim {mustBeNumeric}
        pumpF {mustBeNumeric}
        pumpF2 {mustBeNumeric}
        Ip {mustBeNumeric}
        I0 {mustBeNumeric}
        Idc {mustBeNumeric}
        Istar {mustBeNumeric}
        len {mustBeNumeric}
        modes {mustBeNumeric}
        betanl {mustBeNumeric}
    end 

    methods
        function getK = k(self, freq)
            getK = interp1(self.fsim,self.ksim,freq);
        end
        
        function getG = g(self, freq)
            getG = interp1(self.fsim,self.gsim,freq);
        end

        function self = addMode(self,mode,varargin)
            
            if length(varargin) == 1
                order = varargin{1};
            elseif self.Idc == 0
                order = 4;
            else
                order = 5;
            end

            % Look for instances of pumps
            if isempty(regexpi(mode,'p', 'once'))
                p = 0;
            else
                if isempty(regexpi(mode,'\dp', 'once'))
                    if  isempty(regexpi(mode,'-p', 'once'))
                        p = 1;
                    else
                        p = -1;
                    end
                else
                    if  isempty(regexpi(mode,'-\dp', 'once'))
                        p = str2num(mode(regexpi(mode,'\dp', 'once')));
                    else
                        p = -str2num(mode(regexpi(mode,'\dp', 'once')));
                    end
                end
            end
               
            
            % Look for instances of signals
            if isempty(regexpi(mode,'s', 'once'))
                s = 0;
            else
                if isempty(regexpi(mode,'\ds', 'once'))
                    if  isempty(regexpi(mode,'-s', 'once'))
                        s = 1;
                    else
                        s = -1;
                    end
                else
                    if  isempty(regexpi(mode,'-\ds', 'once'))
                        s = str2num(mode(regexpi(mode,'\ds', 'once')));
                    else
                        s = -str2num(mode(regexpi(mode,'\ds', 'once')));
                    end
                end
            end
            
            % Look for instances of idlers
            if isempty(regexpi(mode,'i', 'once'))
                i = 0;
            else
                if isempty(regexpi(mode,'\di', 'once'))
                    if  isempty(regexpi(mode,'-i', 'once'))
                        i = 1;
                    else
                        i = -1;
                    end
                else
                    if  isempty(regexpi(mode,'-\di', 'once'))
                        i = str2num(mode(regexpi(mode,'\di', 'once')));
                    else
                        i = -str2num(mode(regexpi(mode,'\di', 'once')));
                    end
                end
            end
            

            % add only 3WM idlers
            if order == 3
                self.modes = cat(1,self.modes,[p+i s-i]);
            % add only 4WM idlers
            elseif order == 4
                self.modes = cat(1,self.modes,[p+2*i s-i]);
            % add both idlers
            elseif i ~= 0
                    self.modes = cat(1,self.modes,[p+i s-i]);
                    self.modes = cat(1,self.modes,[p+2*i s-i]);                
            else
                self.modes = cat(1,self.modes,[p s]);
            end
            
        end
    end
end