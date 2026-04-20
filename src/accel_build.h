#ifndef ACCEL_BUILD_H
#define ACCEL_BUILD_H

// Entry point for the `bwa-mem2 build-accel` subcommand. Parses args,
// runs the build pipeline, writes a cache file, returns 0 on success.
int accel_build_main(int argc, char **argv);

#endif // ACCEL_BUILD_H
