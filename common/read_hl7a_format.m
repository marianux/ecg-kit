%% Reads ECG recording in HL7a format
% Reads ECG recordings in Hl7a format from the Chinese database (CCDD). Implements the
% standard "HL7 aECG Implementation Guide - March 21, 2005" available in: 
% 
% https://www.hl7.org/documentcenter/public_temp_75706E59-1C23-BA17-0C25A0CC0545890C/wg/rcrim/annecg/aECG%20Implementation%20Guide%202005-03-21%20final%203.pdf
% 
% Arguments:
%   + filename: recording to be read.
%   + start_sample: (opt) start sample to read. Default 1.
%   + end_sample: (opt) end sample to read. Default min(All recording, ECG block of 200 Mbytes)
% 
% Output:
%   + ECG: the ECG block
%   + heasig: header with the ECG properties. 
%   + ann: annotations for the ECG recordings.
% 
% Limits:
% This routine is limited to read blocks smaller than 200 Mbytes for
% performance reasons. You can disable this limit by doing:
% MaxIOread = Inf; %megabytes
% 
% See also read_ishne_ann, read_ishne_header, read_ECG, ECGwrapper
% 
% Author: Mariano Llamedo Soria
% <matlab:web('mailto:llamedom@electron.frba.utn.edu.ar','-browser') (email)> 
% Version: 0.1 beta
% Birthdate: 02/05/2016
% Last update: 02/05/2016
% Copyright 2008-2016
% 
function [ ECG, heasig, ann, last_sample ] = read_hl7a_format(filename, start_sample, end_sample)


%No leer bloques mas grandes de 200 megabytes
MaxIOread = 200; %megabytes

if( nargin < 2 || isempty( start_sample ) )
    start_sample = 1;
else
    start_sample = max(1,start_sample);
end

ann = [];
heasig = [];
ECG = [];
last_sample = [];

if( nargout > 1 )
    bHeaderRequired = true;
else
    bHeaderRequired = false;
end

if( nargout > 2 )
    bAnnRequired = true;
else
    bAnnRequired = false;
end

try
   xDoc = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

% get series in the file
allSeries = xDoc.getElementsByTagName('series');

if(allSeries.getLength == 0)
    error('read_hl7a_format:NoSeries', 'No series found in %s.\n', filename);
elseif(allSeries.getLength > 1)
    warning('read_hl7a_format:MoreSeries', 'More than one serie in %s. Reading only the first one.\n', filename);
end

allSeries = allSeries.item(0);
etime = allSeries.getElementsByTagName('effectiveTime');
etime = etime.item(0);

loww = etime.getElementsByTagName('low');
loww = loww.item(0);
la = loww.getAttributes;
la = la.item(0);
low_val = char(la.getValue);

loww = etime.getElementsByTagName('high');
loww = loww.item(0);
la = loww.getAttributes;
la = la.item(0);
high_val = char(la.getValue);

sequenceSet = allSeries.getElementsByTagName('sequenceSet');
sequenceSet = sequenceSet.item(0);

allComponents = sequenceSet.getElementsByTagName('component');

for ii = 0:(allComponents.getLength-1)

    thisComp = allComponents.item(ii);
    
    allVals = thisComp.getElementsByTagName('value');
    
    thisVal = allVals.item(0);
    
    thisVal_att = thisVal.getAttributes;
    
    for jj = 0:(thisVal_att.getLength-1)
       
        this_att = thisVal_att.item(jj);

        if( strcmpi(this_att.getValue, 'SLIST_PQ') )
        % ECG leads
            
        elseif( strcmpi(this_att.getValue, 'GLIST_PQ') )
        % time sequence -> sampling rate
            
        end
        
    end
    
end

component = sequenceSet.('component');


sequenceSet

%Note that the item list index is zero-based.
for i=0:allListItems.getLength-1
    thisListItem = allListItems.item(i);
    childNode = thisListItem.getFirstChild;

    while ~isempty(childNode)
        %Filter out text, comments, and processing instructions.
        if childNode.getNodeType == childNode.ELEMENT_NODE
            %Assume that each element has a single org.w3c.dom.Text child
            childText = char(childNode.getFirstChild.getData);
            switch char(childNode.getTagName)
                case 'label' ; itemFound = strcmp(childText,infoLabel);
                case 'callback' ; infoCbk = childText;
            end
        end
        childNode = childNode.getNextSibling;
    end
    if itemFound break; else infoCbk = ''; end
end
disp(sprintf('Item "%s" has a callback of "%s".',infoLabel,infoCbk))


if( nargin < 3 || isempty( end_sample ) )
    %Intento la lectura total por defecto
    samples2read = heasig.nsamp - (start_sample-1);
else
    samples2read = min(heasig.nsamp, end_sample) - (start_sample-1);
end

% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch
   error('Unable to parse XML file %s.',filename);
end

function disp_elem_values( element )

    for ii = 0:(element.getLength-1)

        thisComp = element.item(ii);

    end



% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end
