function graphic_handle = patch_checked(S)
% graphic_handle = patch_checked(S);
% wraps around call to 'patch(S)'.
% checks for empty arguments.
%
% Copyright (C) 2004-2008  Dima Sorkin
% This is free software; and you are welcome to redistribute it
% under certain conditions; see the source for copying conditions.
% There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR
% A PARTICULAR PURPOSE.

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

if ( ~isfield(S,{'vertices'}) )
    error('patch_checked: inp arg does not have a field "vertices"');
end

if ( isempty(S.vertices) )
    warning('PATCH_CHECKED:NO_VERTICES',...
            'patch_checked: inp arg have an empty field "vertices"');
end

error_condition = ( ~isfield(S,{'vertices'}) ) || ( isempty(S.vertices) ) ;

if ( error_condition )
    graphic_handle = [];
    return;
end

graphic_handle = patch(S);
% will not work in 'Octave' this way, only in 'Matlab'.

%%% Octave: handle = patch( % 'CDataMapping', 'scaled',
%%% Octave:                 % 'EdgeColor', 'flat',
%%% Octave:                 % 'FaceVertexCData', S.FaceVertexCData,
%%% Octave:                 'faces', S.faces,
%%% Octave:                 'vertices', S.vertices( : , 1:2 ) );

