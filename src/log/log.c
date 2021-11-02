/**
 * @file log/logger.c
 *
 * @brief Logging library
 *
 * This library can be used to produce logs during simulation runs.
 *
 * SPDX-FileCopyrightText: 2008-2021 HPDCS Group <rootsim@googlegroups.com>
 * SPDX-License-Identifier: GPL-3.0-only
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <arch/io.h>
#include <log/log.h>
#include <core/core.h>

/// The file to write logging information to
static FILE *logfile = NULL;

/// The textual representations and the color codes of the logging levels
static const struct {
	const char *name;
	const char *color;
} levels[] = {
	[LOG_TRACE] = {.name = "TRACE", .color = "\x1b[94m"},
	[LOG_DEBUG] = {.name = "DEBUG", .color = "\x1b[36m"},
	[LOG_INFO] = {.name = "INFO", .color = "\x1b[32m"},
	[LOG_WARN] = {.name = "WARN", .color = "\x1b[33m"},
	[LOG_ERROR] = {.name = "ERROR", .color = "\x1b[31m"},
	[LOG_FATAL] = {.name = "FATAL", .color = "\x1b[35m"}
};

/**
 * @brief Logs a message. Don't use this function directly, rely on the logger() macro instead.
 *
 * @param level the importance level of the message to logger
 * @param file the source file from which the log function was invoked
 * @param line the line number where this function is being called
 * @param fmt a printf-style format string for the message to logger
 * @param ... the list of arguments to fill in the format string @a fmt
 */
void vlogger(unsigned level, char *file, unsigned line, const char *fmt, ...)
{
	va_list args;
	char time_string[IO_TIME_BUFFER_LEN];

	// Check if this invocation should be skipped
	if(level >= global_config.log_level)
		return;

	io_local_time_get(time_string);

	fprintf(
		logfile,
		"%s %s%-5s\x1b[0m \x1b[90m%s:%u:\x1b[0m ",
		time_string,
		levels[level].color,
		levels[level].name,
		file,
		line
	);

	va_start(args, fmt);
	vfprintf(logfile, fmt, args);
	va_end(args);
	fprintf(logfile, "\n");
	fflush(logfile);
}

/**
 * @brief Prints a fancy ROOT-Sim logo
 */
static void print_logo(void)
{
	fprintf(stderr, "\x1b[94m   __ \x1b[90m __   _______   \x1b[94m  _ \x1b[90m       \n");
	fprintf(stderr, "\x1b[94m  /__)\x1b[90m/  ) /  ) /  __ \x1b[94m ( `\x1b[90m . ___ \n");
	fprintf(stderr, "\x1b[94m / \\ \x1b[90m(__/ (__/ (      \x1b[94m._)\x1b[90m / / / )\n");
	fprintf(stderr, "\x1b[0m\n");
}

/**
 * @brief Initialize the logging subsystem
 *
 * @param file The FILE to write logging information to. If set to NULL, it is defaulted to stdout.
 */
void log_init(FILE *file)
{
	logfile = file;
	if(logfile == NULL)
		logfile = stdout;

	if(global_config.log_level != LOG_SILENT)
		print_logo();
}
