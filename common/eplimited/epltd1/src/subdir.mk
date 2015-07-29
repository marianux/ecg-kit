################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/analbeat.c \
../src/bdac.c \
../src/classify.c \
../src/eplimited.c \
../src/match.c \
../src/noisechk.c \
../src/postclas.c \
../src/qrsdet.c \
../src/qrsdet2.c \
../src/qrsfilt.c \
../src/rythmchk.c 

OBJS += \
./src/analbeat.o \
./src/bdac.o \
./src/classify.o \
./src/eplimited.o \
./src/match.o \
./src/noisechk.o \
./src/postclas.o \
./src/qrsdet.o \
./src/qrsdet2.o \
./src/qrsfilt.o \
./src/rythmchk.o 

C_DEPS += \
./src/analbeat.d \
./src/bdac.d \
./src/classify.d \
./src/eplimited.d \
./src/match.d \
./src/noisechk.d \
./src/postclas.d \
./src/qrsdet.d \
./src/qrsdet2.d \
./src/qrsfilt.d \
./src/rythmchk.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<" -I./../../wfdb/windows-amd64/lib/
	@echo 'Finished building: $<'
	@echo ' '


