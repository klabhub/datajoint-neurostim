classdef badBy < dynamicprops
    % A class for managing flags and associated parameters.
    % This class allows dynamic creation of 'badBy...' properties, which
    % are integer vectors, via direct assignment or the setProperty method.
    % It includes a 'parameters' struct and a dependent property 'all'
    % that concatenates unique values from all 'badBy...' data properties.
    % For each 'badBy...' property, a corresponding dependent 'n_badBy...'
    % property is created to count its elements. A 'BadByDataNames' property
    % lists the names of all 'badBy...' data properties.

    properties
        parameters = struct(); % A struct to hold various parameters
    end

    properties (Dependent)
        all % Concatenated unique values of all 'badBy...' data properties
        categories % Cell array of 'badBy...' data property names
    end

    methods
        function disp(self)
            fprintf('BadBy (total: %d unique channels)\n',numel(self.all))
            for c=self.categories'
                fprintf('\t %s: %d [%s ] \n',c,self.("n_" + c),strjoin(string(self.(c)),','))
            end
        end
        function obj = badBy(initialParams)
            % badBy Construct an instance of this class
            if nargin > 0
                if isstruct(initialParams)
                    obj.parameters = initialParams;
                else
                    error('badBy:Constructor:InvalidInput', 'Initial parameters must be a struct.');
                end
            end
        end

        function set.parameters(obj, value)
            if ~isstruct(value)
                error('badBy:setparameters:InvalidType', 'The ''parameters'' property must be a struct.');
            end
            obj.parameters = value;
        end

        function names = get.categories(obj)
            % get.BadByDataNames Getter for the list of 'badBy...' data property names.
            % These are properties that are Dynamic, not Dependent, and start with 'badBy'.
            props = string(properties(obj));
            names = string(props(startsWith(props, 'badBy')));
        end

        function value = get.all(obj)
            % get.all Getter method for the dependent property 'all'.
            % It collects all values from properties listed in BadByDataNames,
            % concatenates them, and returns the unique sorted values.
            value = [];
            dataPropNames = obj.categories; % Use the new dependent property

            for i = 1:length(dataPropNames)
                propName = dataPropNames{i};
                currentBadByData = obj.(propName);
                if ~isempty(currentBadByData) && isnumeric(currentBadByData)
                    value = [value; currentBadByData(:)]; % Ensure column vector
                end
            end
            value = unique(value);
        end

        function setProperty(obj, propName, value)
            % setProperty Sets or creates a 'badBy...' property and its
            % corresponding 'n_badBy...' count property.
            % If propName does not start with 'badBy', it will be
            % transformed: e.g., 'issue' becomes 'badByIssue'.

            if ~ischar(propName) && ~isstring(propName)
                error('badBy:setProperty:InvalidPropertyNameType', 'Property name must be a character array or string.');
            end
            originalPropName = char(propName);
            
            if strcmpi(originalPropName, 'parameters') || strcmpi(originalPropName, 'all') || ...
               startsWith(lower(originalPropName), 'n_badby') || strcmpi(originalPropName, 'badByDataNames')
                error('badBy:setProperty:ReservedName', ...
                      'Cannot use setProperty for reserved names like ''parameters'', ''all'', ''badByDataNames'', or ''n_badBy...''. Assign to ''parameters'' directly if needed.');
            end

            finalPropName = originalPropName;
            if ~startsWith(finalPropName, 'badBy')
                if isempty(finalPropName)
                    error('badBy:setProperty:EmptyPropertyName', 'Property name cannot be empty.');
                end
                finalPropName(1) = upper(finalPropName(1));
                finalPropName = ['badBy', finalPropName];
            end

            if ~isvector(value) || ~isnumeric(value) || ~all(round(value) == value)
                error('badBy:setProperty:InvalidPropertyValue', ...
                      ['Value for ''', finalPropName, ''' (from original input ''', originalPropName, ''') must be an integer vector.']);
            end

            if ~isprop(obj, finalPropName)
                p_data = obj.addprop(finalPropName);
                % p_data.SetObservable = true; % Optional
                % p_data.GetObservable = true; % Optional
            end
            obj.(finalPropName) = int32(value(:));

            countPropertyName = ['n_', finalPropName];
            if ~isprop(obj, countPropertyName)
                p_count = obj.addprop(countPropertyName);
                p_count.Dependent = true;
                p_count.GetMethod = @(~) numel(obj.(finalPropName));
                p_count.SetAccess = 'private';
            else
                metaPropN = findprop(obj, countPropertyName);
                if metaPropN.Dependent
                    metaPropN.GetMethod = @(~) numel(obj.(finalPropName));
                end
            end
        end

        function obj = subsasgn(obj, s, val)
            if isscalar(s) && strcmp(s(1).type, '.')
                propName = s(1).subs;

                if strcmp(propName, 'parameters')
                    obj.parameters = val;
                elseif strcmp(propName, 'all') || startsWith(propName, 'n_badBy') || strcmp(propName, 'badByDataNames')
                    error('MATLAB:class:SetProhibited', ...
                          ['Setting the ''', propName, ''' property of the ''', ...
                           class(obj), ''' class is not allowed.']);
                else
                    obj.setProperty(propName, val);
                end
            else
                obj = builtin('subsasgn', obj, s, val);
            end
        end

        function flg = flag(obj, flg)

            if isscalar(flg)
                flg = repmat("",flg,1);
            end

            cats = obj.categories;

            n_cat = numel(cats);

            for ii = 1:n_cat

                catN = cats(ii);
                flg(obj.(catN)) = catN;

            end

        end
    end

    methods (Static)
        function obj = loadobj(s)
            if isstruct(s)
                obj = badBy();
                fields = fieldnames(s);
                for i = 1:length(fields)
                    fieldName = fields{i};
                    if strcmp(fieldName, 'parameters')
                        obj.parameters = s.parameters;
                    elseif startsWith(fieldName, 'badBy') && ~startsWith(fieldName, 'n_badBy')
                        obj.setProperty(fieldName, s.(fieldName));
                    end
                end
            else
                obj = s;
                mc = metaclass(obj);
                propsList = mc.PropertyList;
                for i = 1:length(propsList)
                    prop = propsList(i);
                    if prop.Dynamic && ~prop.Dependent && startsWith(prop.Name, 'badBy')
                        dataPropName = prop.Name;
                        countPropName = ['n_', dataPropName];
                        if isprop(obj, countPropName)
                            metaCountProp = findprop(obj, countPropName);
                            if metaCountProp.Dynamic && metaCountProp.Dependent
                                metaCountProp.GetMethod = @(~) numel(obj.(dataPropName));
                            else
                                warning('badBy:loadobj:IncorrectCountPropertyState', ...
                                    'Property %s is not correctly configured as a dependent count property.', countPropName);
                            end
                        else
                             warning('badBy:loadobj:MissingCountProperty', ...
                                 'Recreating missing count property: %s', countPropName);
                             p_count = obj.addprop(countPropName);
                             p_count.Dependent = true;
                             p_count.GetMethod = @(~) numel(obj.(dataPropName));
                             p_count.SetAccess = 'private';
                        end
                    end
                end
            end
        end
    end
end