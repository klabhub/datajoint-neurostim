function ip = myip(suffix)
% Function to return the IP address of the computer that runs this function.
% INPUT
% suffix - If provided, the function will reutrn the IP on the network adapter that
% has the specified DNS Suffix.
% OUPUT
% ip  - IP Address as a string
%
% Example
% ip = myip; % Returns all IP addresses on the current machine.
% ip = myip('rad.rutgers.edu');  Returns the ip address that is on the
%                               rad.rutgers.edu network.
%
% BK - Sept 2023

arguments
    suffix = '';
end

if ~ispc;error('This will only work on a Windows machine.');end 

[~,cmdout] = system('ipconfig /all');
if ~isempty(suffix)
    cmdout = extractAfter(cmdout,suffix);
end

ipMatch = regexp(cmdout,'IPv4 Address[\.\s]+:\s+(?<ip>[\d\.]+)','names');

if isempty(ipMatch)    
    ip = "";
else
    ip = string({ipMatch.ip});
end
if ~isempty(suffix)
    ip = ip(1);        
end


