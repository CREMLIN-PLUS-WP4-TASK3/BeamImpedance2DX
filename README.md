# Description

This is a measurement software for particle acceleration structures using bead pulling method.
It utilizes [*PyVISA*](https://pyvisa.readthedocs.io/en/latest/) device communication library as an intermediate device communication layer in order to support multiple platforms.
Device specific code is isolated in driver objects; adding new device driver is a matter of implementing a set of object methods.

This program includes drivers for motors:

- Fake motor (placeholder driver)
- Randomly failing software motor (for program debugging)
- [UIRobot motor driver](http://www.uirobot.com/html/products/944.html)

This program includes drivers for vector analyzers:

- Fake analyzer (placeholder driver)
- Randomly failing software analyzer (for program debugging)
- [Agilent8753et](https://www.keysight.com/en/pd-1000002292%3Aepsg%3Apro-pn-8753ES/)


# Installation

## Get python and pip

Download and install python and pip.

You may use a package manager for this. For example, in Debian installation command looks like this:

`sudo apt-get install python3 python3-pip`

If you don't have a package manager, download installer from the [python web site](https://www.python.org/downloads/).

You will need to install either *NI* drivers (for Windows installations) or *Linux-GPIB* drivers.

## Install *NI* drivers

If you need to use *NI* drivers, use [this manual](https://pyvisa.readthedocs.io/en/latest/faq/getting_nivisa.html#faq-getting-nivisa) for installation.

## Install *Linux-GPIB* drivers

In order to use [*Linux-GPIB*](https://linux-gpib.sourceforge.io/) drivers, one must compile them and install python bindings.

*Note. Drivers are build for current *Linux* kernel. So, after kernel update one needs to recompile GPIB drivers.*

### Getting sources

In order to proceed, you will need to install additional programs:
```bash
apt install -y git git-svn
```

Getting the latest sources:
```bash
git svn clone -r HEAD https://svn.code.sf.net/p/linux-gpib/code/trunk linux-gpib-code
```

Now you have the latest sources of *Linux-GPIB*:
```bash
cd linux-gpib-code
```

### Compiling drivers

For driver compilation you will need Linux sources:
```bash
sudo apt install linux-source
```

```bash
cd linux-gpib-kernel
```

Issue following command to start compilation:

```bash
make
```

And then:

```bash
make install
```

That's it for driver compilation. Return one folder back.

```bash
cd ..
```

### Compiling userspace utilities

```bash
cd linux-gpib-user
```

Since we use source from project repo, we need to prepare [autotools](https://www.gnu.org/software/automake/) installation scripts:

```bash
apt install autotools-dev
./bootstrap
./configure --sysconfdir=/etc
make
make install
```

### Installing Python3 bindings

Go to python bindings folder:
```bash
cd language/python
```

Creating and installing Python3 module:
```bash
python3 setup.py sdist
pip3 install dist/*.tar.gz
```

### Configuring GPIB device

Download GPIB board firmware form <https://linux-gpib.sourceforge.io/firmware/>.

Later in this document `agilent_82350a` device will be used.
Firmware archive is assumed to be located in `~/Downloads/` folder.

Configure your board in `/etc/gpib.conf`. Consult [configuration manual](https://linux-gpib.sourceforge.io/doc_html/configuration-gpib-conf.html).
Configuration will look something like this:
```ini
interface {
        minor = 0       /* board index, minor = 0 uses /dev/gpib0, minor = 1 uses /dev/gpib1, etc. */
        board_type = "agilent_82350b"   /* type of interface board being used */
        name = "violet" /* optional name, allows you to get a board descriptor using ibfind() */
        pad = 0 /* primary address of interface             */
        sad = 0 /* secondary address of interface           */
        timeout = T3s   /* timeout for commands */

        eos = 0x0a      /* EOS Byte, 0xa is newline and 0xd is carriage return
*/
        set-reos = yes  /* Terminate read if EOS */
        set-bin = no    /* Compare EOS 8-bit */
        set-xeos = no   /* Assert EOI whenever EOS byte is sent */
        set-eot = yes   /* Assert EOI with last byte on writes */

/* settings for boards that lack plug-n-play capability */
        base = 0        /* Base io ADDRESS                  */
        irq  = 0        /* Interrupt request level */
        dma  = 0        /* DMA channel (zero disables)      */
}

device {
        minor = 0       /* minor number for interface board this device is connected to */
        name = "analyzer"       /* device mnemonic */
        pad = 7 /* The Primary Address */
        sad = 0 /* Secondary Address */
        timeout = T3s

        eos = 0xa       /* EOS Byte */
        set-reos = no /* Terminate read if EOS */
        set-bin = no /* Compare EOS 8-bit */
}
```

Upload board firmware, if necessary:
```bash
sudo gpib_config --init-data ~/Downloads/gpib_firmware-2008-08-10/hp_82350a/agilent_82350a.bin
```

### Configure udev rules

Some boards require firmware upload on every start.
Also, superuser privileges are required for communication with device.
In order to remove these problems, lets add udev rule for device.

Firstly, lets create new group, which will be able to use GPIB devices:
```bash
addgroup gpib
```

Copy your device's firmware to `/usr/local/sbin`.
Command will look something like this:
```bash
cp Downloads/gpib_firmware-2008-08-10/hp_82350a/agilent_82350a.bin /usr/local/sbin/
```
Then, create a file `/usr/local/sbin/load_agilent` with contents:
```bash
#!/bin/bash

gpib_config --init-data /usr/local/sbin/agilent_82350a.bin
```

And allow its execution:
```bash
chmod 755 /usr/local/sbin/load_agilent
```

Add the following line to the `/etc/udev/rules.d/99-gpib.rules` file:
```udev
KERNEL=="gpib0", SUBSYSTEM=="gpib_common", GROUP="gpib", MODE="0660", RUN+="/usr/local/sbin/load_agilent"
```

Now device will be available for users in group `gpib`, its firmware will be automatically loaded.

## Install *BeadPull*

From the terminal issue the command

`pip3 install beadpull`

for system-wide installation, or

`pip3 install --user beadpull`

for the local installation.

### Add user to device communication groups [Optional]

*Gnu/Linux* distributions restrict access to peripheral devices.

To grant access to serial devices to user, issue command:
```
adduser user dialout
```
USB devices:
```
adduser user plugdev
```
GPIB devices (if configured):
```
adduser user gpib
```

# Run

In order to run the program, from the terminal issue the command

`beadpull`

If it's not found, you have to modify the `PATH` environment variable to include the python scripts folder ([unix guide](https://stackoverflow.com/a/3402176/5745120),[windows guide](https://superuser.com/a/143121/736971)) or supply the full path to the program.

# Other

## Update *beadpull*

To update program, from the terminal issue the command

`pip3 install --upgrade beadpull`

for system-wide installation, or

`pip3 install --upgrade --user beadpull`

for the local installation.

## Reset *beadpull* settings to defaults

Program saves its settings under
`APPDATA/beadupll` on Windows and `$XDG_CONFIG_HOME/beadpull` or `$HOME/.config/beadpull` on other systems.
In order to reset settings, simply delete saved settings folder.

## Issue submission

In case of a problem with the program, create an issue in the [issue tracker](https://gitlab.com/matsievskiysv/beadpull/issues).
