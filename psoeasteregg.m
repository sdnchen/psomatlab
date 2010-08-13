function state = psoeasteregg(options,state,flag)

if strcmp(flag,'done') && ~strcmp(options.Display,'off')
    sounddir = [getenv('ProgramFiles') '\Starcraft\Sound\Zerg\Advisor'] ;
    if isdir(sounddir) && exist([sounddir '\ZAdUpd02.wav'],'file')
        [y, Fs, nbits] = wavread([sounddir '\ZAdUpd02.wav']) ;
        obj = audioplayer(y, Fs, nbits);
        playblocking(obj)
    end
end