# xmldig2csv Converter

`xmldig2csv` is a command-line tool that converts `.XMLdig` files—containing base64-encoded digital channel data captured with a **Teledyne LeCroy Serial Analyzer (SDA)** and an **HDA125 digitizer**—into human-readable `.csv` files.

The `.XMLdig` file format is briefly documented in the *MAUI Oscilloscopes Remote Control and Automation Manual*, page 200[^1]:

> While there is no corresponding template file for digital waveforms, they are organized similarly to analog waveforms, with a text-based description of the header, followed by a Base64 encoding of the waveform binary data. 
In Base64 encoding, 4 ASCII characters (6 bits per character) are used to represent 3 bytes, which, in the case of the digital waveforms, correspond to the value of 3 bits in the waveform. For example, the initial 4 bytes “AAEA” in the XML data field map to the 3 bits “010”. The following image shows how a small section of the waveform data would translate.


The output format is identical to the CSV format generated by Teledyne's [WaveStudio](https://www.teledynelecroy.com/support/softwaredownload/wavestudio.aspx) software, including a timestamped header and logic values.

## Example Output

```
Time,CS,CA0,CA1,CA2,CA3,CA4,CA5,CA6,CA7,CA8,CA9
-5.024226e-06,1,1,1,1,0,0,1,1,1,1
-5.024146e-06,1,1,1,1,0,0,1,1,1,1
...
```

The labels used for the column headers correspond to the names defined in the digital bus of the TELEDYNE HDA125 digitizer.

## Usage

Run the tool by passing the `.XMLdig` file as argument:

```
./xmldig2csv ~/my-trace.XMLdig
```

This will produce `my-trace.csv` in the same directory. The program is silent and produces no stdout output on success.

## Build Instructions

The project uses CMake and depends on:

- **C++17**
- **LibXML2** (for parsing XML)
- **OpenMP** (optional, for parallelization)

### On Ubuntu/Debian

Install dependencies:

```
sudo apt update
sudo apt install build-essential cmake libxml2-dev libxml++2.6-dev
```

Build the project:

```
git clone https://github.com/mcsee-artifacts/xmldig2csv.git
cd xmldig2csv
mkdir build && cd build
cmake ..
make
```

The resulting binary `xmldig2csv` will be located in the `build/` directory.

> [!IMPORTANT]
> - On Linux, both `libxml2` and `libxml++2.6` are required.
> - If OpenMP is not available, the program falls back to single-threaded processing.
> - The build assumes a Unix-like environment.

## Testing

We provide sample files for testing a working setup in the `example-data/` directory.

[^1]: [MAUI Oscilloscopes Remote Control and Automation Manual (page 200)](https://cdn.teledynelecroy.com/files/manuals/maui-remote-control-and-automation-manual.pdf)
