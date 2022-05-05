# Digital Speech Decoder - Florida Man Edition 
This version of DSD is a flavor blend of [szechyjs](https://github.com/szechyjs/dsd "szechyjs") RTL branch and some of my own additions, along with portions of DMR and NXDN code from the [LouisErigHerve](https://github.com/LouisErigHerve/dsd "LouisErigHerve") branch as well. NXDN voice decoding is currently working a lot better, thanks to the latter, although I have yet to explore the expanded NXDN or DMR decoding he has laid out. That is a goal. I have also implemented a few more RTL options, including rtl device gain, PPM error, device index selection, squelch, VFO bandwidth, and a UDP remote that works like the old rtl_udp fork, although its currently limited to changing frequency and squelch. The goal is to integrate this project into [EDACS-FM](https://github.com/lwvmobile/edacs-fm "EDACS-FM") but I also want it to be its own standalone project. 

## 2022.05.05 Update ## 
I have successfully added a DMR Stereo method for listening to voice audio in both TDMA channels/slots simultaneously. This method will also allow for data decoding in the opposite slot if only one voice call is active, allowing the user not to miss any useful information in the second slot while the previous slot is in use. DMR Stereo also has vastly improved handling of MS/Simplex voice decoding and (limited/hit or miss) MS data decoding. To call this method, see the example below. The old method of decoding DMR is also included and is stil the default for the time being, but the DMR Stereo method is working well, users are encouraged to try both methods and find the one that is suitable for them. New commands include:
```
-T Enable DMR 'TDMA' Stereo
   Please note, for the time being, any MBE file saving, WAV file saving, or MBE payloads are disabled while using DMR    Stereo. These features are still present in the default DMR handling and other voice decoder options and will be re-    integrated into DMR Stereo in the future. DMR Stereo will also need to be enabled to handle DMR MS/Simplex voice, as    it was removed from the default method due to poor performance.
   
-F Enable Passive Frame Sync
   Use Passive Frame Sync if voice in both slots becomes choppy/skips
   Using Passive Sync may cause wonky audio though, depends on quality of signal.
  
```

![DSD-FME](https://github.com/lwvmobile/dsd-fme/blob/pulseaudio/dsd-fme.png)

![DSD-FME](https://github.com/lwvmobile/dsd-fme/blob/pulseaudio/dsd-fme2.png)

## Example Usage - New DMR Stereo, Ncurses Terminal, Pulse Input/Output, and Log Console to file
`./dsd -fr -T -Z -N 2> voice.log`

in second terminal, same folder, run:

`tail -n 40 -f voice.log`

## Example Usage - RTL
`dsd -fi -i rtl -c 154.9875M -P -2 -D 1 -G 36 -L 70 -U 6021 -Y 12`

```
-i rtl to use rtl_fm (default is -i pulse for pulse audio)

-o pulse to set output to pulse audio (is default anyways)

-c Set frequency

-P set PPM error

-D set device index number

-G set device gain (0-49) (default = 0 Auto Gain)

-L set rtl squelch to 25

-V set RTL sample 'volume' multiplier

-U set UDP port for rtl_fm remote control

-Y 12 set rtl VFO bandwidth in kHz, (default = 48)(6, 8, 12, 16, 24, 48)

```

### Note

`dsd` and `./dsd` may need to be used in your usage, depending on if you are running it from the build folder, or if you have run make install. You will need to tweak the examples to your particular usage


## Example Pulse Audio Input and Pulse Audio Output, Autodetect Frame Type
`dsd` or `./dsd`
```
Yes, it really is this simple now

-fa Auto-detect frame type, is default

sans NXDN, dPMR, or Provoice, need to specify -fp, -fd, -fi, or -fn respectively
pulse audio is set as default input and output methods

```
## Example STDIN UDP from GQRX or SDR++, output to Pulse Audio, and save wav files
`socat stdio udp-listen:7355 | dsd -fi -i - -o pulse -w nxdn.wav`

## Example Usage - Extra Information for Academic Study and Logging
`dsd -d ./MBE/ -Z -pu 2>> voice.log`

in second terminal, same folder, run:

`tail -n 40 -f voice.log`

## Roadmap
The Current list of objectives include:

1. ~~Random Tinkering~~ More Random Tinkering

2. ~~Implemented Pulse Audio~~ Remove PortAudio, ~~improved Pulse Audio for stereo output and channel/slot separation~~

3. ~~Improve NXDN and DMR support~~ Continue to improve ~~NXDN and DMR~~ all data and voice decoding. 

4. ~~More Concise Printouts - Ncurses~~

5. ~~Make simple to use installer script.~~ 

## How to clone, check out, and build this branch

### Ubuntu 22.04/20.04/LM20/Debian Bullseye or Newer:

Using the included install.sh should make for a simple and painless clone, build, and install on newer Debian or Ubuntu based system. Simply acquire or copy the script, and run it. Update: Ubuntu 22.04 has been tested working with the installer script and use. 

```
wget https://raw.githubusercontent.com/lwvmobile/dsd-fme/pulseaudio/install.sh
chmod +x install.sh
./install.sh
```

### Ubuntu 18.04/LM19/Buster Note:
The above install.sh should now function on older system types. You can elect to manually follow the steps down below. Do NOT Manually build and install ITPP 4.3.1 on older systems, it is currently not wanting to build on Ubuntu 18.04 and Linux Mint 19. Install it from the repository instead.

## Manual Install

First, install dependency packages. This guide will assume you are using Ubuntu 20.04 based distros. Check your package manager for equivalent packages if different. PortAudio is not currently used in this build, and is disabled in CMakeLists.txt, you can re-enable it if you wish, but it isn't recommended unless you have a very specific reason to do so. Some of these dependencies are not currently be used, but may be used in future builds.


```
sudo apt update
sudo apt install libpulse-dev pavucontrol libsndfile1-dev libfftw3-dev liblapack-dev socat libusb-1.0-0-dev libncurses5 libncurses5-dev rtl-sdr librtlsdr-dev libusb-1.0-0-dev cmake git wget make build-essential libitpp-dev libncursesw5-dev
```


### Build and Install ITPP - ONLY IF NOT IN REPO!!

```
wget -O itpp-latest.tar.bz2 http://sourceforge.net/projects/itpp/files/latest/download?source=files
tar xjf itpp*
#if you can't cd into this folder, double check folder name first
cd itpp-4.3.1 
mkdir build
cd build
cmake ..
make -j `nproc`
sudo make install
sudo ldconfig
cd ..
cd ..
```

MBELib is considered a requirement in this build. You must read this notice prior to continuing. [MBElib Patent Notice](https://github.com/lwvmobile/mbelib#readme "MBElib Patent Notice") This version of MBELib is 1.3.1 and prints to STDERR, using the stock 1.3.0 MBElib may cause problems with print alignments if paired with this version of DSD-FME, or cause issues with printing over the screen in future builds.

```
git clone https://github.com/lwvmobile/mbelib
cd mbelib
mkdir build
cd build
cmake ..
make -j `nproc`
sudo make install
sudo ldconfig
cd ..
cd ..
```

Finish by running these steps to clone and build DSD-FME w/ pulseaudio support.

```
git clone https://github.com/lwvmobile/dsd-fme
cd dsd-fme
git branch -a
git checkout remotes/origin/pulseaudio
git checkout -b pulseaudio
git branch -a #double check to see if you are on pulseaudio branch
mkdir build
cd build
cmake ..
make -j `nproc`
##only run make install if you don't have another version already installed##
sudo make install
sudo ldconfig

```
Optional 'Virtual Sinks' for routing audio from SDR++ or GQRX, Media Players, etc. into DSD-FME

You may wish to direct sound into DSD-FME via Virtual Sinks. You may set up a Virtual Sink or two on your machine for routing audio in and out of applications to other applications using the following command, and opening up pavucontrol "PulseAudio Volume Control" in the menu to change line out of application to virtual sink, and line in of DSD-FME to monitor of virtual sink. This command will not persist past a reboot, so you will need to invoke them each time you reboot, or search for how to add this to your conf files for persistency if desired.

```
pacmd load-module module-null-sink sink_name=virtual_sink  sink_properties=device.description=Virtual_Sink
pacmd load-module module-null-sink sink_name=virtual_sink2  sink_properties=device.description=Virtual_Sink2
```

Already have this branch, and just want to pull the latest build?

```
##Open your clone folder##
git pull https://github.com/lwvmobile/dsd-fme pulseaudio
##cd into your build folder##
cd build
##cmake usually isn't necesary, but could be if I update the cmakelist.txt
cmake ..
make -j `nproc`
sudo make install
sudo ldconfig
```

## License
Copyright (C) 2010 DSD Author
GPG Key ID: 0x3F1D7FD0 (74EF 430D F7F2 0A48 FCE6  F630 FAA2 635D 3F1D 7FD0)

    Permission to use, copy, modify, and/or distribute this software for any
    purpose with or without fee is hereby granted, provided that the above
    copyright notice and this permission notice appear in all copies.

    THE SOFTWARE IS PROVIDED "AS IS" AND ISC DISCLAIMS ALL WARRANTIES WITH
    REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
    AND FITNESS.  IN NO EVENT SHALL ISC BE LIABLE FOR ANY SPECIAL, DIRECT,
    INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
    LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
    OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
    PERFORMANCE OF THIS SOFTWARE.
    
   
