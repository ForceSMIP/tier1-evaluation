function [field1,field2,field3,field4,field5,field6,field7,field8,field9,field10] = get_avg_field_nd(files, varnam_cell)

% GET_AVG_FIELD_ND  Reads fields from a set of netCDF files and averages them.
%  
%   FIELD = GET_AVG_FIELD(FILES, VARNAM) reads the fields with
%   variable name VARNAM (can be a cell array of variable names)
%   from the netCDF files given in the cell FILES. The average 
%   of these fields is returned as FIELD1-FIELD10. Output is only 
%   given up to the nuber of variables specified in VARNAM. 
%   Maximum number of variables is 10

%   Output format is in the format the variables are stored in netcdf

field1 = 0;
field2 = 0;
field3 = 0;
field4 = 0;
field5 = 0;
field6 = 0;
field7 = 0;
field8 = 0;
field9 = 0;
field10 = 0;

if iscell(varnam_cell)
for i = 1:length(varnam_cell) % loop over variables in VARNAM
    varnam = varnam_cell{i};
 
  if isstr(files)    % convert to cell array
    files    = cellstr(files);
  end
  
  % number of files used in averaging
  nfiles     = length(files);
  
  %find alias for pre-fram simulations.
 % varname = alias(files, varnam);
  varname = varnam;
  
  % read fields from files and average them
  sum_field  = 0;
  sum_nan  = 0;
  for ifile = 1:nfiles
    ncid     = netcdf.open(char(files{ifile}), 'NC_NOWRITE');
    field    = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(varname))));
    
     % replace missing values by NaN 
%     miss_val = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,char(varnam)),0);
%     if ~isempty(miss_val)
%       imissing	= find( abs(field - miss_val) <= sqrt(eps));
%              field(imissing) = NaN;
%     end

    field(abs(field)>1e30) = nan;
    sum_nan = sum_nan + isnan(field);
    field(isnan(field)) = 0;
    sum_field = sum_field + field;  
    netcdf.close(ncid);  
  end
  frac_nans = sum_nan ./nfiles;
  field = sum_field ./ (1-frac_nans) ./ nfiles;
  
switch i
    case 1
        field1 = field;
    case 2
        field2 = field;
    case 3
        field3 = field;
    case 4
        field4 = field;
    case 5
        field5 = field;
    case 6
        field6 = field;
    case 7
        field7 = field;
    case 8
        field8 = field;
    case 9
        field9 = field;
    case 10
        field10 = field;
    otherwise
        disp('get_avg_field gets up to 10 variables, all subsequent ignored')
        break
end
  
end
else
  varnam = varnam_cell;
  
  if isstr(files)    % convert to cell array
    files    = cellstr(files);
  end
  
  % number of files used in averaging
  nfiles     = length(files);
  
  %find alias for pre-fram simulations.
 % varname = alias(files, varnam);
  varname = varnam;
  
  % read fields from files and average them
  sum_field  = 0;
  sum_nan  = 0;
  for ifile = 1:nfiles
    ncid     = netcdf.open(char(files{ifile}), 'NC_NOWRITE');
    field    = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(varname))));
    
     % replace missing values by NaN 
%     miss_val = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,char(varnam)),0);
%     if ~isempty(miss_val)
%       imissing	= find( abs(field - miss_val) <= sqrt(eps));
%              field(imissing) = NaN;
%     end

    field(abs(field)>1e30) = nan;
    sum_nan = sum_nan + isnan(field);
    field(isnan(field)) = 0;
    sum_field = sum_field + field;  
    netcdf.close(ncid);  
  end
  frac_nans = sum_nan ./nfiles;
  field = sum_field ./ (1-frac_nans) ./ nfiles;
  field1 = field;
end
      