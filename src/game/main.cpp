#include "main.h"

constexpr Interface::FullscreenMode fullscreen_mode = Interface::borderless_fullscreen;

Interface::Window window(std::string(window_name), screen_size * 2, Interface::windowed, adjust_(Interface::WindowSettings{}, min_size = screen_size));
static Graphics::DummyVertexArray dummy_vao = nullptr;

static Audio::Context context = nullptr;
Audio::SourceManager audio_manager;

Graphics::ShaderConfig shader_config = Graphics::ShaderConfig::Core();
// static Interface::ImGuiController gui_controller(Poly::derived<Interface::ImGuiController::GraphicsBackend_Modern>, adjust_(Interface::ImGuiController::Config{}, shader_header = shader_config.common_header, store_state_in_file = ""));

namespace Fonts
{
    namespace Files
    {
        Graphics::FontFile main("assets/CatIV15.ttf", 15);
    }

    Graphics::Font main;
}

Graphics::TextureAtlas texture_atlas = []{
    Graphics::TextureAtlas ret(ivec2(2048), "assets/_images", "assets/atlas.png", "assets/atlas.refl");
    auto font_region = ret.Get("font_storage.png");

    Unicode::CharSet glyph_ranges;
    glyph_ranges.Add(Unicode::Ranges::Basic_Latin);

    Graphics::MakeFontAtlas(ret.GetImage(), font_region.pos, font_region.size, {
        {Fonts::main, Fonts::Files::main, glyph_ranges, Graphics::FontFile::hinting_mode_light},
    });
    return ret;
}();
Graphics::Texture texture_main = Graphics::Texture(nullptr).Wrap(Graphics::clamp).Interpolation(Graphics::nearest).SetData(texture_atlas.GetImage());

AdaptiveViewport adaptive_viewport(shader_config, screen_size);
Render r = adjust_(Render(0x2000, shader_config), SetTexture(texture_main), SetMatrix(adaptive_viewport.GetDetails().MatrixCentered()));

Input::Mouse mouse;

Random<> randomize(std::random_device{}());

struct ProgramState : Program::DefaultBasicState
{
    GameUtils::State::StateManager state_manager;
    GameUtils::FpsCounter fps_counter;

    void Resize()
    {
        adaptive_viewport.Update();
        mouse.SetMatrix(adaptive_viewport.GetDetails().MouseMatrixCentered());
    }

    Metronome metronome = Metronome(60);

    Metronome *GetTickMetronome() override
    {
        return &metronome;
    }

    int GetFpsCap() override
    {
        return 60 * NeedFpsCap();
    }

    void EndFrame() override
    {
        fps_counter.Update();
        window.SetTitle(STR((window_name), " TPS:", (fps_counter.Tps()), " FPS:", (fps_counter.Fps())));
    }

    void Tick() override
    {
        window.ProcessEvents();
        // window.ProcessEvents({gui_controller.EventHook()});

        if (window.Resized())
        {
            Resize();
            Graphics::Viewport(window.Size());
        }
        if (window.ExitRequested())
            Program::Exit();

        if ((Input::Button(Input::l_alt).down() || Input::Button(Input::r_alt).down()) && Input::Button(Input::enter).pressed())
            window.SetMode(window.Mode() == Interface::windowed ? fullscreen_mode : Interface::windowed);

        // gui_controller.PreTick();
        state_manager.Tick();
        audio_manager.Tick();
    }

    void Render() override
    {
        // gui_controller.PreRender();
        adaptive_viewport.BeginFrame();
        state_manager.Render();
        adaptive_viewport.FinishFrame();
        Graphics::CheckErrors();
        // gui_controller.PostRender();

        window.SwapBuffers();
    }


    void Init()
    {
        // ImGui::StyleColorsDark();

        // // Load various small fonts
        // auto monochrome_font_flags = ImGuiFreeTypeBuilderFlags_LightHinting;

        // gui_controller.LoadFont("assets/CatIV15.ttf", 15.0f, adjust(ImFontConfig{}, FontBuilderFlags = monochrome_font_flags));
        // gui_controller.LoadDefaultFont();
        // gui_controller.RenderFontsWithFreetype();

        Audio::LoadMentionedFiles(Audio::LoadFromPrefixWithExt("assets/sounds/"), Audio::mono, Audio::wav);

        Graphics::Blending::Enable();
        Graphics::Blending::FuncNormalPre();

        mouse.HideCursor();

        state_manager.NextState().Set("Menu");
    }
};

int _main_(int, char **)
{
    ProgramState program_state;
    program_state.Init();
    program_state.Resize();
    program_state.RunMainLoop();
    return 0;
}
