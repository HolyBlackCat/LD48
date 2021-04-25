#pragma once

constexpr ivec2 screen_size = ivec2(480, 270);
constexpr std::string_view window_name = "EARTHWYRM";

extern Interface::Window window;
extern Graphics::ShaderConfig shader_config;

extern Audio::SourceManager audio_manager;

namespace Fonts
{
    namespace Files
    {
        extern Graphics::FontFile main;
    }

    extern Graphics::Font main;
}

extern Graphics::TextureAtlas texture_atlas;
extern Graphics::Texture texture_main;

extern AdaptiveViewport adaptive_viewport;
extern Render r;

extern Input::Mouse mouse;

extern Random<> randomize;

extern GameUtils::FpsCounter fps_counter;
