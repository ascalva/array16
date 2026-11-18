cmake_minimum_required(VERSION 3.13)

# pico-sdk 2.x: pass -DPICO_SDK_PATH=... or set the environment variable.
if (NOT PICO_SDK_PATH)
    set(PICO_SDK_PATH $ENV{PICO_SDK_PATH})
endif ()
if (NOT PICO_SDK_PATH)
    message(FATAL_ERROR "Point PICO_SDK_PATH at a pico-sdk 2.x checkout (environment or -DPICO_SDK_PATH=...)")
endif ()
include(${PICO_SDK_PATH}/external/pico_sdk_import.cmake)

project(array16_controller C CXX ASM)
set(CMAKE_C_STANDARD 11)
pico_sdk_init()

add_executable(array16_controller rp2040/main.c core/ctrl.c)
pico_generate_pio_header(array16_controller ${CMAKE_CURRENT_LIST_DIR}/rp2040/sequencer.pio)
target_include_directories(array16_controller PRIVATE core)
target_compile_options(array16_controller PRIVATE -Wall -Wextra)
target_link_libraries(array16_controller pico_stdlib hardware_pio hardware_dma hardware_clocks hardware_adc)

# The controller speaks a binary protocol over USB CDC: no newline translation, no waiting for a host.
target_compile_definitions(array16_controller PRIVATE PICO_STDIO_DEFAULT_CRLF=0 PICO_STDIO_USB_CONNECT_WAIT_TIMEOUT_MS=0)
pico_enable_stdio_usb(array16_controller 1)
pico_enable_stdio_uart(array16_controller 0)
pico_set_program_name(array16_controller "array16 beam controller")
pico_add_extra_outputs(array16_controller)
