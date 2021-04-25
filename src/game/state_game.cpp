#include "game/main.h"
#include "gameutils/tiled_map.h"
#include "signals/connection.h"
#include "utils/json.h"

SIMPLE_STRUCT( Atlas
    DECL(Graphics::TextureAtlas::Region)
        worm, tiles, vignette, panel, button, button_icons, button_subscripts, cursor, mouse_pointer, bug
)
static Atlas atlas;

constexpr int tile_size = 12;

namespace Draw
{
    struct ShiftedRegion
    {
        const Graphics::TextureAtlas::Region &main, *back = nullptr, *front = nullptr;
        bool unclamp = false;

        ShiftedRegion(const Graphics::TextureAtlas::Region &main)
            : main(main), back(&main), front(&main)
        {}

        ShiftedRegion(const Graphics::TextureAtlas::Region &new_main, const Graphics::TextureAtlas::Region *new_back, const Graphics::TextureAtlas::Region *new_front, bool flip = false, bool new_unclamp = false)
            : main(new_main), back(new_back), front(new_front), unclamp(new_unclamp)
        {
            if (flip)
                std::swap(back, front);
        }
    };

    void ShiftedTile(ivec2 pos, ShiftedRegion region, ivec2 shift, float alpha = 1, float beta = 1)
    {
        if (region.unclamp)
            r.iquad(pos + shift, region.main).alpha(alpha).beta(beta);
        else
            r.iquad(pos + clamp_min(shift, 0), region.main.region(clamp_min(-shift, 0), region.main.size - abs(shift))).alpha(alpha).beta(beta);

        if (shift)
        {
            auto shift_reg = (shift > 0).any() ? region.back : region.front;
            if (shift_reg)
                r.iquad(pos + mod_ex(clamp_max(shift, 0), region.main.size), shift_reg->region(mod_ex(-clamp_min(shift, 0), region.main.size), abs(shift) + (1 - abs(sign(shift))) * region.main.size)).alpha(alpha).beta(beta);
        }
    }

    void Tooltip(std::string_view str)
    {
        constexpr int spacing = 4;
        constexpr ivec2 bg_padding(2, 0);
        ivec2 pos(0, screen_size.y / 2 - atlas.panel.size.y - spacing);
        Graphics::Text text(Fonts::main, str);
        ivec2 text_size = text.ComputeStats().size;
        r.iquad(pos - bg_padding - text_size with(x /= 2), text_size + bg_padding * 2).color(fvec3(0)).alpha(0.75);
        r.itext(pos, text).align(fvec2(0,1)).color(fvec3(1));
    }

    void InitAtlas()
    {
        static bool once = true;
        if (once)
        {
            once = false;
            texture_atlas.InitRegions(atlas, ".png");
        }
    }

    void MousePointer()
    {
        if (window.HasMouseFocus())
            r.iquad(mouse.pos(), atlas.mouse_pointer).center();
    }
}

namespace Sounds
{
    #define ADD_SOUND(name_, var_) \
        Audio::Source &name_(fvec2 pos, float vol = 1, float pitch = 0) \
        { \
            return audio_manager.Add(Audio::File(#name_##_c))->pos(pos).volume(vol).pitch(pow(2, pitch + float(-var_ <= randomize.real() <= var_))).play(); \
        } \
        Audio::Source &name_(float vol = 1, float pitch = 0) \
        { \
            return audio_manager.Add(Audio::File(#name_##_c))->relative().volume(vol).pitch(pow(2, pitch + float(-var_ <= randomize.real() <= var_))).play(); \
        } \

    ADD_SOUND(dirt_breaks, 0.5)
    ADD_SOUND(dirt_placed, 0.2)
    ADD_SOUND(landing, 0.2)
    ADD_SOUND(click, 0.03)
    ADD_SOUND(player_dies, 0.1)
    ADD_SOUND(bug_dies, 0.3)
    ADD_SOUND(explosion, 0.2)
    ADD_SOUND(reverse, 0.2)
    ADD_SOUND(bait, 0.2)
    ADD_SOUND(next_level, 0.2)

    #undef ADD_SOUND
}

namespace Music
{
    Audio::Buffer buffer;
    Audio::Source source;

    bool enabled = true;

    void StartOnce()
    {
        static bool once = true;
        if (once)
        {
            once = false;
            buffer = Audio::Sound(Audio::ogg, Audio::stereo, "assets/theme.ogg");
            source = buffer;
            source.volume(2.5);
            source.loop();
            if (enabled)
                source.play();
        }
    }

    void Toggle()
    {
        enabled = !enabled;
        if (source)
        {
            if (enabled)
                source.play();
            else
                source.pause();
        }
    }

    void ToggleIfHotkeyPressed()
    {
        if (Input::Button(Input::m).pressed())
            Toggle();
    }
}

namespace DebugInfo
{
    static bool enabled =
        #ifdef NDEBUG
        false;
        #else
        true;
        #endif

    void Tick()
    {
        if (Input::Button(Input::f1).pressed())
            enabled = !enabled;
    }

    void Render()
    {
        if (!enabled)
            return;
        r.itext(-screen_size/2, Graphics::Text(Fonts::main, FMT("TPS: {}\nFPS: {}\n[F1] to hide this.", fps_counter.Tps(), fps_counter.Fps()))).align(fvec2(-1)).color(fvec3(1)).alpha(0.5);
    }
}

namespace Gui
{
    class ButtonList;

    struct Button
    {
        Sig::Connection<Button, ButtonList> con;

        int icon_index = 0;
        std::string tooltip = 0;
        int hover_timer = 0;

        bool force_blocked = false;
        std::optional<int> subscript;
        virtual bool IsBlocked() const {return force_blocked || (subscript && *subscript <= 0);}

        Button() {}
        Button(int icon_index, std::string tooltip) : icon_index(icon_index), tooltip(std::move(tooltip)) {}

        ivec2 pos{};

        static int Size()
        {
            return atlas.button.size.y;
        }
        static int SideDecoW()
        {
            static int ret = atlas.button.size.x % atlas.button.size.y;
            return ret;
        }

        virtual int GetIconIndex() const
        {
            return icon_index;
        }

        virtual bool IsHovered() const
        {
            if (IsBlocked())
                return false;
            return (mouse.pos() >= pos).all() && (mouse.pos() < pos + Size()).all();
        }

        virtual bool IsClicked() const
        {
            if (IsBlocked())
                return false;
            return mouse.left.released() && IsHovered();
        }

        virtual void Tick()
        {
            if (IsClicked())
                Click();

            if (IsHovered())
                clamp_var_max(++hover_timer, 1000);
            else
                hover_timer = 0;
        }

        virtual void Click()
        {
            Sounds::click(mouse.pos());
        }

        virtual int GetBackgroundVariant() const
        {
            return IsHovered() ? (mouse.left.down() ? 2 : 1) : 0;
        }

        virtual void Draw() const
        {
            int bg_variant = GetBackgroundVariant();
            r.iquad(pos, atlas.button.region(ivec2(SideDecoW() + Button::Size() * bg_variant, 0), ivec2(Button::Size())));
            r.iquad(pos, atlas.button_icons.region(ivec2(Button::Size() * GetIconIndex(), 0), ivec2(Button::Size()))).alpha(bg_variant == 2 ? 0.75 : 1);
            if (subscript)
            {
                r.iquad(pos + Size() - ivec2(atlas.button_subscripts.size.y), atlas.button_subscripts.region(ivec2(atlas.button_subscripts.size.y * clamp(*subscript, 0, 9), 0), ivec2(atlas.button_subscripts.size.y)));
                if (IsBlocked())
                    r.iquad(pos, ivec2(Size())).tex(atlas.button.pos + atlas.button.size - 0.5f, fvec2(0)).alpha(0.5);
            }
            if (hover_timer > 30)
                Draw::Tooltip(tooltip);
        }
    };

    struct Checkbox : Button
    {
        int enabled_icon_index = 0;
        bool enabled = false;
        Checkbox() {}
        Checkbox(int icon_index, int enabled_icon_index, bool enabled, std::string tooltip) : Button(icon_index, std::move(tooltip)), enabled_icon_index(enabled_icon_index), enabled(enabled) {}

        int GetIconIndex() const override
        {
            return enabled ? enabled_icon_index : icon_index;
        }

        void Click() override
        {
            Button::Click();
            enabled = !enabled;
        }
    };

    struct StickyButton : Button
    {
        using Button::Button;

        bool was_clicked = false;

        bool IsHovered() const override
        {
            return !was_clicked && Button::IsHovered();
        }
        bool IsClicked() const override
        {
            return !was_clicked && Button::IsClicked();
        }

        void Click() override
        {
            Button::Click();
            was_clicked = true;
        }
    };

    struct RadioButton;
    struct RadioButtonList
    {
        Sig::ConnectionList<RadioButtonList, RadioButton> con;
    };

    struct RadioButton : Button
    {
        Sig::Connection<RadioButton, RadioButtonList> other_buttons;

        bool active = false;

        RadioButton(int icon_index, RadioButtonList &others, std::string tooltip)
            : Button(icon_index, std::move(tooltip))
        {
            Sig::Bind<&RadioButton::other_buttons, &RadioButtonList::con>(*this, others);
        }

        int GetBackgroundVariant() const override
        {
            return Button::GetBackgroundVariant() + 3 * active;
        }

        void Tick() override
        {
            if (IsBlocked())
                active = false;
            Button::Tick();
        }

        void Click() override
        {
            bool activate = !active;

            const auto &list = other_buttons->con;
            for (size_t i = 0; i < list.Size(); i++)
                list[i]->active = false;

            active = activate;
            Button::Click();
        }
    };

    class ButtonList
    {
      public:
        int x = 0;
        float align = 0;

        Sig::ConnectionList<ButtonList, Button> con;

        ButtonList() {}
        ButtonList(int x, float align, const std::vector<Button *> &buttons)
            : x(x), align(align)
        {
            for (Button *button : buttons)
                Sig::Bind<&Button::con, &ButtonList::con>(*button, *this);

        }

        size_t NumButtons() const
        {
            return con.Size();
        }

        ivec2 CalcFirstButtonPos() const
        {
            return ivec2(x - NumButtons() * Button::Size() * (align * 0.5 + 0.5), screen_size.y / 2 - Button::Size());
        }

        void Tick()
        {
            ivec2 first_pos = CalcFirstButtonPos();
            for (size_t i = 0; i < NumButtons(); i++)
            {
                con[i]->pos = first_pos with(x += Button::Size() * i);
                con[i]->Tick();
            }
        }

        void Draw() const
        {
            ivec2 pos = CalcFirstButtonPos();

            // Side decorations.
            r.iquad(pos with(x -= Button::SideDecoW()), atlas.button.region(ivec2(0), ivec2(Button::SideDecoW(), atlas.button.size.y)));
            r.iquad(pos with(x += NumButtons() * Button::Size()), atlas.button.region(ivec2(0), ivec2(Button::SideDecoW(), atlas.button.size.y))).flip_x();

            for (size_t i = 0; i < NumButtons(); i++)
                con[i]->Draw();
        }
    };
}

struct Particle
{
    fvec2 pos{}, vel{}, acc{};
    float damp = 0;
    int time = 0;
    int max_time = 30;

    SIMPLE_STRUCT_WITHOUT_NAMES( Params
        DECL(float INIT=4) size
        DECL(fvec3 INIT=fvec3(1)) color
        DECL(float INIT=1) alpha
        DECL(float INIT=1) beta
        VERBATIM Params() = default;
    )
    struct TimedParams : Params
    {
        int t = 0;
        TimedParams() {}
        TimedParams(const Params &other) : Params(other) {}
    };

    Params p;
    std::optional<TimedParams> p1, p2;

    struct ParamRange
    {
        const Params *a = nullptr;
        const Params *b = nullptr;
        float ta = 0;
        float tb = 0;

        auto Calc(int time, auto &&func) const
        {
            if (!b)
                return func(a);
            return mix((time - ta) / float(tb - ta), func(a), func(b));
        }
    };

    ParamRange GetParamRange() const
    {
        ParamRange ret;
        if (!p1 && !p2)
        {
            ret.a = &p;
            return ret;
        }
        if (p1.has_value() != p2.has_value())
        {
            ret.a = &p;
            ret.ta = 0;
            ret.b = p1 ? &*p1 : &*p2;
            ret.tb = p1 ? p1->t : p2->t;
            return ret;
        }
        const TimedParams *pa = &*p1, *pb = &*p2;
        if (pa->t > pb->t)
            std::swap(pa, pb);
        if (time < pa->t)
        {
            ret.a = &p;
            ret.ta = 0;
            ret.b = pa;
            ret.tb = pa->t;
            return ret;
        }
        else
        {
            ret.a = pa;
            ret.ta = pa->t;
            ret.b = pb;
            ret.tb = pb->t;
            return ret;
        }
    }

    void Draw(ivec2 viewport_center, ivec2 viewport_size) const
    {
        fvec2 rel_pos = pos - viewport_center;
        if ((abs(rel_pos) > viewport_size / 2 + 16).any())
            return;

        auto range = GetParamRange();

        Params params;
        Meta::cexpr_for<Refl::Class::member_count<Params>>([&](auto index)
        {
            constexpr auto i = index.value;
            Refl::Class::Member<i>(params) = range.Calc(time, [](const Params *p){return Refl::Class::Member<i>(*p);});
        });

        r.fquad(rel_pos, fvec2(params.size)).center().color(params.color).alpha(params.alpha).beta(params.beta);
    }
};

struct ParticleController
{
    std::deque<Particle> particles;

    void Tick()
    {
        std::erase_if(particles, [](const Particle &pa)
        {
            return pa.time >= pa.max_time;
        });

        for (Particle &pa : particles)
        {
            pa.pos += pa.vel;
            pa.vel += pa.acc;
            pa.vel *= (1 - pa.damp);
            pa.time++;
        }
    }

    void Render(ivec2 viewport_center, ivec2 viewport_size) const
    {
        for (const Particle &pa : particles)
            pa.Draw(viewport_center, viewport_size);
    }

    static Particle ParticleDirt(fvec2 pos, fvec2 vel, float min_size, float max_size, int max_time)
    {
        float p = 0 <= randomize.real() <= 1;
        return adjust(Particle{}
            , pos = pos
            , vel = vel
            , acc = fvec2(0,0.05)
            , max_time = max_time
            , damp = mix(p, 0.01, 0.03)
            , p.size = mix(p, max_size, min_size)
            , p.color = (std::array{fvec3(139, 90, 60) / 255, fvec3(101, 62, 41) / 255, fvec3(159, 111, 56) / 255}[randomize.integer() < 3] * float(0.9 <= randomize.real() <= 1.1))
            , p1 = _object_.p
            , p1->t = _object_.max_time
            , p1->size = 0
        );
    }

    static Particle ParticleBlood(fvec2 pos, fvec2 vel, float min_size, float max_size, int max_time)
    {
        float p = 0 <= randomize.real() <= 1;
        return adjust(Particle{}
            , pos = pos
            , vel = vel
            , acc = fvec2(0,0.05)
            , max_time = max_time
            , damp = mix(p, 0.01, 0.03)
            , p.size = mix(p, max_size, min_size)
            , p.color = fvec3(0.4 <= randomize.real() <= 1, 0, 0)
            , p1 = _object_.p
            , p1->t = _object_.max_time
            , p1->size = 0
        );
    }

    static Particle ParticleFlame(fvec2 pos, fvec2 vel, float min_size, float max_size, int max_time)
    {
        float p = 0 <= randomize.real() <= 1;
        float c = -0.5 <= randomize.real() <= 1;
        return adjust(Particle{}
            , pos = pos
            , vel = vel
            , acc = fvec2(0,-0.03)
            , max_time = max_time
            , damp = mix(p, 0.03, 0.05)
            , p.size = mix(p, max_size, min_size)
            , p.color = clamp(fvec3(c + 1, c, 0))
            , p1 = _object_.p
            , p1->t = _object_.max_time
            , p1->size = 0
        );
    }

    static Particle ParticleBait(fvec2 pos, fvec2 vel, float min_size, float max_size, int max_time)
    {
        float p = 0 <= randomize.real() <= 1;
        float c = -0.3 <= randomize.real() <= 0.3;
        return adjust(Particle{}
            , pos = pos
            , vel = vel
            , acc = fvec2(0,0.02)
            , max_time = max_time
            , damp = mix(p, 0.06, 0.12)
            , p.size = mix(p, max_size, min_size)
            , p.color = clamp(fvec3(1 - abs(c), c, -c))
            , p1 = _object_.p
            , p1->t = _object_.max_time
            , p1->beta = 0
            , p1->size = 0
        );
    }


    void EffectDirtDestroyed(fvec2 center_pos, fvec2 area, int n = 35)
    {
        while (n-- > 0)
        {
            particles.push_back(ParticleDirt(
                center_pos + fvec2(-area.x/2 <= randomize.real() <= area.x/2, -area.y/2 <= randomize.real() <= area.y/2),
                fvec2::dir(randomize.angle(), 0 <= randomize.real() <= 1),
                2, 8, 30 <= randomize.integer() <= 60)
            );
        }
    }

    void EffectDirtMinor(fvec2 center_pos, fvec2 area, int n)
    {
        while (n-- > 0)
        {
            particles.push_back(ParticleDirt(
                center_pos + fvec2(-area.x/2 <= randomize.real() <= area.x/2, -area.y/2 <= randomize.real() <= area.y/2),
                fvec2::dir(randomize.angle(), 0 <= randomize.real() <= 1),
                2, 6, 20 <= randomize.integer() <= 40)
            );
        }
    }

    void EffectBloodExplosion(fvec2 center_pos, fvec2 area, int n = 5)
    {
        while (n-- > 0)
        {
            particles.push_back(ParticleBlood(
                center_pos + fvec2(-area.x/2 <= randomize.real() <= area.x/2, -area.y/2 <= randomize.real() <= area.y/2),
                fvec2::dir(randomize.angle(), 0 <= randomize.real() <= 1.3),
                2, 12, 40 <= randomize.integer() <= 70)
            );
        }
    }

    void EffectBloodExplosionMinor(fvec2 center_pos, fvec2 area, int n = 5)
    {
        while (n-- > 0)
        {
            particles.push_back(ParticleBlood(
                center_pos + fvec2(-area.x/2 <= randomize.real() <= area.x/2, -area.y/2 <= randomize.real() <= area.y/2),
                fvec2::dir(randomize.angle(), 0 <= randomize.real() <= 2),
                2, 8, 15 <= randomize.integer() <= 45)
            );
        }
    }

    void EffectExplosion(fvec2 center_pos, fvec2 area, int n = 20)
    {
        while (n-- > 0)
        {
            particles.push_back(ParticleFlame(
                center_pos + fvec2(-area.x/2 <= randomize.real() <= area.x/2, -area.y/2 <= randomize.real() <= area.y/2),
                fvec2::dir(randomize.angle(), 0 <= randomize.real() <= 1.25),
                2, 12, 15 <= randomize.integer() <= 30)
            );
        }
    }

    void EffectBait(fvec2 center_pos, fvec2 area, int n = 20)
    {
        while (n-- > 0)
        {
            particles.push_back(ParticleBait(
                center_pos + fvec2(-area.x/2 <= randomize.real() <= area.x/2, -area.y/2 <= randomize.real() <= area.y/2),
                fvec2::dir(randomize.angle(), 0.5 <= randomize.real() <= 2),
                2, 6, 30 <= randomize.integer() <= 60)
            );
        }
    }

    void EffectBaitMinor(fvec2 center_pos, fvec2 area)
    {
        EffectBait(center_pos, area, 10);
    }
};

namespace MapPoints
{
    struct BugSpawn
    {
        ivec2 pos{}, dir{};
    };

    ENUM( Direction,
        (right)(down)(left)(up)
    )

    ivec2 DirFromEnum(Direction value)
    {
        return ivec2::dir4(int(value));
    }

    struct State
    {
        bool worm_placed = false;
        ivec2 worm_starting_pos{};
        ivec2 worm_starting_dir{};
        int worm_starting_len = 4;

        std::string level_name;

        bool exit_placed = false;
        ivec2 exit_pos{};
        ivec2 exit_dir{};

        std::vector<BugSpawn> bug_spawns;

        int num_dirt = 0;
        int num_bait = 0;
        int num_bombs = 0;
        int num_reverses = 0;
    };

    STRUCT( Base POLYMORPHIC )
    {
        UNNAMED_MEMBERS()
        virtual void Process(State &state, ivec2 tile_pos) = 0;
    };

    STRUCT( Start EXTENDS Base )
    {
        MEMBERS(
            DECL(std::string ATTR Refl::Optional) name
            DECL(int) len
            DECL(Direction) dir
            DECL(int INIT=0 ATTR Refl::Optional) dirt
            DECL(int INIT=0 ATTR Refl::Optional) bait
            DECL(int INIT=0 ATTR Refl::Optional) bombs
            DECL(int INIT=0 ATTR Refl::Optional) rev
        )

        void Process(State &state, ivec2 tile_pos) override
        {
            if (state.worm_placed)
                Program::Error("More than one worm.");
            state.worm_placed = true;

            state.level_name = name;

            state.worm_starting_pos = tile_pos;
            state.worm_starting_dir = DirFromEnum(dir);
            state.worm_starting_len = len;
            state.num_dirt = dirt;
            state.num_bait = bait;
            state.num_bombs = bombs;
            state.num_reverses = rev;
        }
    };

    STRUCT( Exit EXTENDS Base )
    {
        MEMBERS(
            DECL(Direction) dir
        )

        void Process(State &state, ivec2 tile_pos) override
        {
            if (state.exit_placed)
                Program::Error("More than one exit.");
            state.exit_placed = true;

            state.exit_pos = tile_pos;
            state.exit_dir = DirFromEnum(dir);
        }
    };

    STRUCT( Bug EXTENDS Base )
    {
        MEMBERS(
            DECL(Direction) dir
        )

        void Process(State &state, ivec2 tile_pos) override
        {
            state.bug_spawns.push_back({tile_pos, DirFromEnum(dir)});
        }
    };
}

class Map
{
  public:
    // Tile types.
    enum class Tile
    {
        air,
        dirt,
        pipe,
        wall,
        spike,
        bait,
        _count,
    };

    enum class TileRenderer {simple, fancy, pipe};

    struct TileInfo
    {
        TileRenderer renderer{};
        int tile_index = -1;
        bool breakable = true;
        int path_priority = 0; // Worm prefers tiles with a higher value. 0 means it will refuse to go there.
        bool kills = false;
        bool bombable = false;
    };

    // Tile properties table.
    inline static const TileInfo tile_info[] =
    {
        /* air   */ {.renderer = TileRenderer::simple, .tile_index = -1, .breakable = false, .path_priority = 99, .kills = false, .bombable = false},
        /* dirt  */ {.renderer = TileRenderer::fancy , .tile_index =  0, .breakable = true , .path_priority =  9, .kills = false, .bombable = true },
        /* pipe  */ {.renderer = TileRenderer::pipe  , .tile_index =  1, .breakable = false, .path_priority =  0, .kills = false, .bombable = false},
        /* wall  */ {.renderer = TileRenderer::pipe  , .tile_index =  2, .breakable = false, .path_priority =  0, .kills = false, .bombable = true },
        /* spike */ {.renderer = TileRenderer::pipe  , .tile_index =  3, .breakable = false, .path_priority =  9, .kills = true , .bombable = true },
        /* bait  */ {.renderer = TileRenderer::simple, .tile_index =  0, .breakable = true , .path_priority = 99, .kills = false, .bombable = true },
    };
    static_assert(std::size(tile_info) == size_t(Tile::_count));

    static const TileInfo &Info(Tile tile)
    {
        ASSERT(int(tile) < int(Tile::_count), "Tile index is out of range.");
        return tile_info[int(tile)];
    }

    // Tile merging rules.
    static bool ShouldMergeWith(Tile a, Tile b)
    {
        if (a == Tile::spike)
            return b != Tile::spike && b != Tile::air;

        if (a == b)
            return true;

        if (a == Tile::dirt && b == Tile::wall)
            return true;
        return false;
    }

    struct TileData
    {
        Tile type = Tile::air;
    };

    uint64_t global_time = 0;

    Array2D<TileData> tiles;
    int level_index = -1;

    Array2D<uint8_t> noise;

    MapPoints::State data;

    inline static std::string map_prefix = "assets/levels/";

    static int GetLevelCount()
    {
        static int ret =
        []{
            int ret = 0;
            for (auto elem : std::filesystem::directory_iterator(map_prefix))
            {
                if (elem.path().extension() != ".json")
                    continue;
                std::string name = elem.path().stem().string();
                if (std::find_if_not(name.begin(), name.end(), Stream::Char::IsDigit{}) != name.end())
                    continue; // If the name contains anything but digits, skip this file.
                int index = Refl::FromString<unsigned short>(name);
                clamp_var_min(ret, index);
            }
            return ret;
        }();
        return ret;
    }

    Map() {}
    Map(int new_level_index)
        : level_index(new_level_index)
    {
        try
        {
            if (level_index >= GetLevelCount())
                Program::Error("Map index is out of range.");

            Json json(Stream::ReadOnlyData::file(FMT("{}{}.json", map_prefix, new_level_index+1)).string(), 32);

            // Parse tiles.
            Array2D<int> layer_mid = Tiled::LoadTileLayer(Tiled::FindLayer(json.GetView(), "mid"));
            tiles = decltype(tiles)(layer_mid.size());
            noise = decltype(noise)(layer_mid.size());
            for (index_vec2 pos : vector_range(tiles.size()))
            {
                int index = layer_mid.safe_nonthrowing_at(pos);
                if (index < 0 || index >= int(Tile::_count))
                    Program::Error(FMT("Tile index {} at position {} is out of range.", index, pos));
                TileData new_data;
                new_data.type = Tile(index);
                tiles.safe_nonthrowing_at(pos) = new_data;
                noise.safe_nonthrowing_at(pos) = randomize.integer<decltype(noise)::type>();
            }

            // Parse points.
            Tiled::PointLayer points = Tiled::LoadPointLayer(Tiled::FindLayer(json.GetView(), "obj"));
            for (const auto &[name, pos] : points.points)
            {
                ivec2 tile_pos = div_ex(iround(pos), tile_size);
                Refl::FromString<Refl::PolyStorage<MapPoints::Base>>(name)->Process(data, tile_pos);
            }
        }
        catch (std::exception &e)
        {
            Program::Error(STR("While loading map #", (new_level_index+1), ":\n", (e.what())));
        }
    }

    MAYBE_CONST(
        CV TileData &At(ivec2 pos) CV
        {
            return tiles.clamped_at(pos);
        }
    )

    const TileInfo &InfoAt(ivec2 pos) const
    {
        return Info(At(pos).type);
    }

    decltype(noise)::type Noise(ivec2 pos) const
    {
        return noise.safe_nonthrowing_at(mod_ex(pos, noise.size()));
    }

    bool TilePosIsAtBorder(ivec2 pos) const
    {
        return (pos - 1 < 0).any() || (pos + 1 >= tiles.size()).any();
    }

    void render(ivec2 viewport_center, ivec2 viewport_size) const
    {
        ivec2 corner_a = div_ex(viewport_center - viewport_size / 2, tile_size);
        ivec2 corner_b = div_ex(viewport_center + viewport_size / 2, tile_size);
        for (ivec2 tile_pos : corner_a <= vector_range <= corner_b)
        {
            ivec2 pixel_pos = tile_pos * tile_size - viewport_center;
            const TileData& data = At(tile_pos);
            Tile type = data.type;
            const TileInfo info = tile_info[int(type)];

            switch (info.renderer)
            {
              case TileRenderer::simple:
                if (info.tile_index < 0)
                    continue;
                r.iquad(pixel_pos, atlas.tiles.region(ivec2(0, info.tile_index) * tile_size, ivec2(tile_size)));
                break;
              case TileRenderer::fancy:
                {
                    if (info.tile_index < 0)
                        continue;

                    bool merge_full = true;
                    for (int i = 0; i < 8; i++)
                    {
                        if (!ShouldMergeWith(type, At(tile_pos + ivec2::dir8(i)).type))
                        {
                            merge_full = false;
                            break;
                        }
                    }

                    Graphics::TextureAtlas::Region base_region = atlas.tiles.region(ivec2(1, 2 * info.tile_index) * tile_size, ivec2(tile_size));

                    if (merge_full)
                    {
                        int variants[] = {0,0,0,0,0,1,1,1,2,3};
                        int var = variants[Noise(tile_pos) % std::size(variants)];
                        base_region.pos += bvec2(var & 1, var & 2) * tile_size;
                        r.iquad(pixel_pos, base_region);
                        break;
                    }

                    auto corner = [&](ivec2 sub)
                    {
                        ivec2 tile_offset = sub * 2 - 1;
                        ivec2 pixel_offset = sub * (tile_size / 2);

                        bool merge_h = ShouldMergeWith(type, At(tile_pos + tile_offset with(y = 0)).type);
                        bool merge_v = ShouldMergeWith(type, At(tile_pos + tile_offset with(x = 0)).type);
                        bool merge_d = ShouldMergeWith(type, At(tile_pos + tile_offset            ).type);

                        Graphics::TextureAtlas::Region region;

                        if (merge_h && merge_v && merge_d)
                        {
                            region = base_region;
                        }
                        else if (merge_h && merge_v)
                        {
                            region = base_region with(pos.x += tile_size * 3);
                        }
                        else if (merge_h)
                        {
                            region = base_region with(pos.x += tile_size * 2);
                        }
                        else if (merge_v)
                        {
                            region = base_region with(pos += tile_size * ivec2(2,1));
                        }
                        else
                        {
                            region = base_region with(pos += tile_size * ivec2(3,1));
                        }

                        r.iquad(pixel_pos + pixel_offset, region.region(pixel_offset, ivec2(tile_size / 2)));
                    };

                    corner(ivec2(0,0));
                    corner(ivec2(0,1));
                    corner(ivec2(1,0));
                    corner(ivec2(1,1));
                }
                break;
              case TileRenderer::pipe:
                {
                    int merge = 0;
                    for (int i = 0; i < 4; i++)
                        merge = merge << 1 | ShouldMergeWith(type, At(tile_pos + ivec2::dir4(i)).type);
                    ivec2 variant(0);
                    if (merge == 0b0010)
                        variant = ivec2(2, 0);
                    else if (merge == 0b0001)
                        variant = ivec2(1, 1);
                    else if (merge == 0b1000)
                        variant = ivec2(2, 1);
                    else if (merge == 0b0100)
                        variant = ivec2(1, 0);
                    else if (merge == 0b1010)
                        variant = ivec2(0, 0);
                    else if (merge == 0b0101)
                        variant = ivec2(0, 1);
                    else if (merge == 0b1100)
                        variant = ivec2(3, 0);
                    else if (merge == 0b0110)
                        variant = ivec2(4, 0);
                    else if (merge == 0b0011)
                        variant = ivec2(4, 1);
                    else if (merge == 0b1001)
                        variant = ivec2(3, 1);
                    else
                        variant = ivec2(5, 0);
                    r.iquad(pixel_pos, atlas.tiles.region((ivec2(1, 2 * info.tile_index) + variant) * tile_size, ivec2(tile_size)));
                }
                break;
            }
        }
    }
};

namespace Progress
{
    constexpr std::string_view file_path = "progress.txt";

    SIMPLE_STRUCT( Data
        DECL(int INIT=0) level // 1-based.
    )

    [[nodiscard]] Data Load()
    {
        try
        {
            return Refl::FromString<Data>(Stream::Input(file_path), adjust(Refl::FromStringOptions{}, ignore_missing_fields = true));
        }
        catch (...)
        {
            return {};
        }
    }

    void Save(const Data &data)
    {
        try
        {
            Stream::Output output = std::string(file_path);
            Refl::ToString(data, output);
            output.Flush();
        }
        catch (...)
        {
            // Do nothing.
        }
    }
}

struct Worm
{
    std::vector<ivec2> segments; // This is sorted tail-to-head.
    int crawl_offset = 0;
    float fall_offset = 0;
    float fall_speed = 0;
    bool dead = false;
    bool need_death_anim = false;
    bool no_sound_on_death_animation = false;
    int death_timer = 0;
    bool out_of_bounds = false;

    ivec2 ChooseDirection(const Map &map) const
    {
        if (segments.size() < 2)
            return ivec2();
        ivec2 head = segments.back();
        ivec2 pre_head = segments[segments.size() - 2];
        ivec2 current_dir = head - pre_head;
        return ChooseDirectionLow(map, head, current_dir, [&](ivec2 tile_pos)
        {
            return std::any_of(segments.begin(), segments.end(), [&](ivec2 seg){return seg == tile_pos;});
        });
    }

    static ivec2 ChooseDirectionLow(const Map &map, ivec2 pos, ivec2 current_dir, std::function<bool(ivec2 tile_pos)> is_tile_banned = nullptr, bool go_back_if_stuck = false, bool prefer_left = false)
    {
        if (!current_dir)
            current_dir = map.data.worm_starting_dir;

        auto GetPathPriority = [&](ivec2 tile_pos) -> int
        {
            int ret = 0;

            if (is_tile_banned && is_tile_banned(tile_pos))
                ret = 0; // Refuse to touch your own body.
            else if (tile_pos == map.data.exit_pos && tile_pos - pos == -map.data.exit_dir)
                ret = 98; // Go the exit if possible.
            else
                ret = map.InfoAt(tile_pos).path_priority;

            return ret;
        };

        auto p_fwd = GetPathPriority(pos + current_dir);
        auto p_left = GetPathPriority(pos + current_dir.rot90(-1));
        auto p_right = GetPathPriority(pos + current_dir.rot90(1));
        auto p_back = go_back_if_stuck ? GetPathPriority(pos - current_dir) : 0;

        // Nowhere to go.
        if (max(p_fwd, p_left, p_right) <= 0)
        {
            if (p_back > 0)
                return -current_dir;
            return ivec2();
        }

        // If the current dir is at least as good as the alternatives, keep it.
        if (p_fwd >= p_left && p_fwd >= p_right)
            return current_dir;

        // Check for bait.
        if (p_left == p_right)
        {
            int min_dist_sqr = -1;
            std::optional<ivec2> bait_dir;
            for (index_vec2 bait_tile : vector_range(map.tiles.size()))
            {
                if (map.At(bait_tile).type == Map::Tile::bait && abs(sign(bait_tile - pos)) != abs(current_dir) )
                {
                    int dist_sqr = (bait_tile - pos).len_sqr();
                    if (min_dist_sqr == -1 || min_dist_sqr > dist_sqr)
                    {
                        min_dist_sqr = dist_sqr;
                        bait_dir = bait_tile - pos;
                    }
                }
            }
            if (bait_dir)
                return sign(abs(sign(current_dir.rot90())) * *bait_dir);
        }

        // If `go_back_if_stuck` is set, and both left and right are equally good, go back.
        if (go_back_if_stuck && p_left == p_right)
            return -current_dir;

        // Pick the best of the two directions.
        if (prefer_left ? p_left >= p_right : p_left > p_right)
            return current_dir.rot90(-1);
        else
            return current_dir.rot90(1);
    }

    void Tick(Map &map, ParticleController &par, bool unpaused)
    {
        // Death conditions.
        if (!dead)
        {
            // Spikes.
            for (ivec2 seg : segments)
            {
                if (map.InfoAt(seg).kills)
                {
                    dead = true;
                    need_death_anim = true;
                    break;
                }
            }
        }

        // Death animation.
        if (need_death_anim)
        {
            need_death_anim = false;
            for (ivec2 seg : segments)
                par.EffectBloodExplosion(seg * tile_size + tile_size/2, fvec2(tile_size));

            if (!no_sound_on_death_animation)
                Sounds::player_dies(segments[segments.size() / 2] * tile_size + tile_size / 2);
        }

        // Stop if dead.
        if (dead)
        {
            clamp_var_max(++death_timer, 1000);
            return;
        }

        // Stop if paused.
        if (!unpaused)
            return;

        bool on_ground = false;
        std::vector<unsigned char> supported_segments(segments.size());
        { // Check if the worm is supported by the ground.
            auto TileIsSolid = [&](ivec2 pos)
            {
                const Map::TileInfo &data = map.InfoAt(pos);
                return data.path_priority < 99 && !data.kills;
            };

            std::optional<ivec2> prev_delta, next_delta;
            for (size_t i = 0; i < segments.size(); i++) LOOP_NAME(ground)
            {
                if (TileIsSolid(segments[i] with(y++)))
                {
                    on_ground = true;
                    supported_segments[i] = true;
                    continue;
                }

                prev_delta = next_delta;
                next_delta = i == segments.size() - 1 ? std::nullopt : std::optional{segments[i+1] - segments[i]};

                for (std::optional<ivec2> *delta_ptr : {&prev_delta, &next_delta})
                {
                    if (!*delta_ptr)
                        continue;
                    for (int s : {-1, 1})
                    {
                        ivec2 offset = (*delta_ptr)->rot90(s);
                        if (offset != ivec2(0,-1) && TileIsSolid(segments[i] + offset))
                        {
                            on_ground = true;
                            supported_segments[i] = true;
                            break;
                        }
                    }
                    if (supported_segments[i])
                        break;
                }
            }
        }

        if (on_ground)
        {
            if (fall_speed > 0.8)
            {
                Sounds::landing(0.2);
                for (size_t i = 0; i < segments.size(); i++)
                {
                    if (supported_segments[i])
                        par.EffectDirtMinor(segments[i] * tile_size + ivec2(tile_size / 2, tile_size), fvec2(tile_size, 0), 4);
                }
            }
            fall_offset = 0;
            fall_speed = 0;
        }
        else
        {
            clamp_var_max(fall_speed += 0.02, 1);
            fall_offset += fall_speed;
            if (fall_offset >= tile_size)
            {
                fall_offset = 0;
                for (ivec2 &seg : segments)
                    seg.y++;
            }
        }

        // Crawl.
        if (on_ground && segments.size() >= 2)
        {
            int crawl_dist = map.global_time % 2;

            // We crawl pixel by pixel, just in case.
            while (crawl_dist-- > 0)
            {
                bool ok = true;
                if (crawl_offset >= tile_size / 2 - 3)
                {
                    // Play the winning sound.
                    if (segments.back() == map.data.exit_pos && segments[segments.size() - 2] != map.data.exit_pos)
                        Sounds::next_level(map.data.exit_pos * tile_size + tile_size/2);

                    bool at_exit = map.data.exit_pos == segments.back();
                    ivec2 chosen_dir = at_exit ? ivec2() : ChooseDirection(map);

                    if (chosen_dir || at_exit)
                    {
                        crawl_offset -= tile_size;
                        segments.front() = segments.back() + chosen_dir;
                        std::rotate(segments.begin(), segments.begin() + 1, segments.end());

                        { // Destroy tiles.
                            if (!map.TilePosIsAtBorder(segments.back()) && map.InfoAt(segments.back()).breakable)
                            {
                                Map::Tile prev_tile = map.At(segments.back()).type;
                                map.At(segments.back()).type = Map::Tile::air;
                                ivec2 pixel_pos = segments.back() * tile_size + tile_size / 2;
                                if (prev_tile == Map::Tile::bait)
                                {
                                    par.EffectBait(pixel_pos, fvec2(tile_size));
                                    Sounds::bait(pixel_pos, 0.5);
                                }
                                else
                                {
                                    par.EffectDirtDestroyed(pixel_pos, fvec2(tile_size));
                                    Sounds::dirt_breaks(pixel_pos, 0.5);
                                }
                            }
                        }

                        // Check if out of bounds.
                        if (!map.tiles.pos_in_range(segments.back()))
                            out_of_bounds = true;
                    }
                    else
                    {
                        ok = false;
                    }
                }
                if (ok)
                    crawl_offset++;
            }
        }
    }

    void Draw(ivec2 camera_pos) const
    {
        float alpha = clamp_min(1 - death_timer / 30.f);
        if (alpha <= 0)
            return;

        std::optional<ivec2> prev_delta, next_delta;
        for (size_t i = 0; i < segments.size(); i++)
        {
            prev_delta = next_delta;
            next_delta = i == segments.size() - 1 ? std::nullopt : std::optional{segments[i+1] - segments[i]};

            Graphics::TextureAtlas::Region region = atlas.worm.region(ivec2(0), ivec2(tile_size));

            ivec2 base_pos = segments[i] * tile_size - camera_pos + ivec2(0, fall_offset);

            if ((prev_delta && !*prev_delta) || (next_delta && !*next_delta))
                continue;

            if (!prev_delta && !next_delta)
            {
                // This shouldn't happen, draw a placeholder.
                r.iquad(base_pos, region).alpha(alpha);
            }
            else if (!prev_delta && next_delta)
            {
                // Tail.
                if (next_delta == ivec2(1,0))
                    region = region.region(ivec2(2,2) * tile_size, ivec2(tile_size));
                else if (next_delta == ivec2(-1,0))
                    region = region.region(ivec2(3,2) * tile_size, ivec2(tile_size));
                else if (next_delta == ivec2(0,-1))
                    region = region.region(ivec2(2,3) * tile_size, ivec2(tile_size));
                else
                    region = region.region(ivec2(3,3) * tile_size, ivec2(tile_size));

                Graphics::TextureAtlas::Region straight_reg;
                if (crawl_offset)
                {
                    if (next_delta->x != 0)
                        straight_reg = atlas.worm.region(ivec2(0,0) * tile_size, ivec2(tile_size));
                    else
                        straight_reg = atlas.worm.region(ivec2(1,0) * tile_size, ivec2(tile_size));
                }

                Draw::ShiftedTile(base_pos, {region, nullptr, &straight_reg, (*next_delta < 0).any(), true}, *next_delta * crawl_offset, alpha);
            }
            else if (prev_delta && !next_delta)
            {
                // Head.
                if (prev_delta == ivec2(-1,0))
                    region = region.region(ivec2(0,1) * tile_size, ivec2(tile_size));
                else if (prev_delta == ivec2(1,0))
                    region = region.region(ivec2(1,1) * tile_size, ivec2(tile_size));
                else if (prev_delta == ivec2(0,1))
                    region = region.region(ivec2(0,2) * tile_size, ivec2(tile_size));
                else
                    region = region.region(ivec2(1,2) * tile_size, ivec2(tile_size));

                Graphics::TextureAtlas::Region straight_reg;
                if (crawl_offset)
                {
                    if (prev_delta->x != 0)
                        straight_reg = atlas.worm.region(ivec2(0,0) * tile_size, ivec2(tile_size));
                    else
                        straight_reg = atlas.worm.region(ivec2(1,0) * tile_size, ivec2(tile_size));
                }

                Draw::ShiftedTile(base_pos, {region, segments.size() > 2 ? &straight_reg : nullptr, nullptr, (*prev_delta < 0).any(), true}, *prev_delta * crawl_offset, alpha);
            }
            else if (*prev_delta == *next_delta)
            {
                // Straight piece.
                if (next_delta->x != 0)
                    region = region.region(ivec2(0,0) * tile_size, ivec2(tile_size));
                else
                    region = region.region(ivec2(1,0) * tile_size, ivec2(tile_size));

                Draw::ShiftedTile(base_pos, {region, i == 1 ? nullptr : &region, i == segments.size() - 2 ? nullptr : &region, (*prev_delta < 0).any() }, *next_delta * crawl_offset, alpha);
            }
            else
            {
                // Corner piece.
                ivec2 delta = prev_delta->rot90() == *next_delta ? *prev_delta : -*next_delta;
                if (delta == ivec2(0,-1))
                    region = region.region(ivec2(2,0) * tile_size, ivec2(tile_size));
                else if (delta == ivec2(1,0))
                    region = region.region(ivec2(3,0) * tile_size, ivec2(tile_size));
                else if (delta == ivec2(-1,0))
                    region = region.region(ivec2(2,1) * tile_size, ivec2(tile_size));
                else
                    region = region.region(ivec2(3,1) * tile_size, ivec2(tile_size));
                r.iquad(base_pos, region).alpha(alpha);
            }
        }
    }
};

struct Bug
{
    ivec2 tile_pos;
    ivec2 dir;

    int crawl = 0;

    bool dead = false;
    int death_timer = 0;

    static constexpr int max_death_timer = 20;

    ivec2 HitboxTilePos() const
    {
        if (crawl < tile_size / 2)
            return tile_pos;
        else
            return tile_pos + dir;
    }

    ivec2 ChooseDirection(const Map &map) const
    {
        return Worm::ChooseDirectionLow(map, tile_pos, dir, [&](ivec2 pos){return map.InfoAt(pos).path_priority < 99 && map.At(pos).type != Map::Tile::spike;}, true, true);
    }

    void Tick(Map &map, ParticleController &par)
    {
        // Crawl.
        if (!dead && map.global_time % 2 == 0)
        {
            if (crawl > 0)
            {
                crawl++;
                if (crawl == tile_size)
                {
                    crawl = 0;
                    tile_pos += dir;
                }
            }
            else
            {
                if (ivec2 next_dir = ChooseDirection(map))
                {
                    dir = next_dir;
                    crawl++;
                }
            }
        }

        // Destroy the bait if touching it.
        if (!dead)
        {
            ivec2 hitbox_pos = HitboxTilePos();
            if (map.At(hitbox_pos).type == Map::Tile::bait)
            {
                map.At(hitbox_pos).type = Map::Tile::air;
                ivec2 pixel_pos = hitbox_pos * tile_size + tile_size / 2;
                par.EffectBait(pixel_pos, fvec2(tile_size));
                Sounds::bait(pixel_pos, 0.5);
            }
        }

        // Die from spikes.
        if (!dead && map.InfoAt(HitboxTilePos()).kills)
        {
            ivec2 pixel_pos = HitboxTilePos() * tile_size + tile_size / 2;
            par.EffectBloodExplosionMinor(pixel_pos, fvec2(tile_size));
            Sounds::bug_dies(pixel_pos, 0.3);
            dead = true;
        }

        if (dead)
            death_timer++;
    }

    void Render(ivec2 camera_pos) const
    {
        constexpr int anim_period = 40;
        constexpr int sprite_size = tile_size * 2;

        int anim_state = window.Ticks() % anim_period < anim_period / 2;
        int anim_dir = dir == ivec2(1,0) ? 0 : dir == ivec2(0,1) ? 1 : dir == ivec2(-1,0) ? 2 : 3;
        float alpha = clamp(1 - death_timer / float(max_death_timer));

        r.iquad(tile_pos * tile_size + tile_size / 2 + dir * crawl - camera_pos, atlas.bug.region(ivec2(anim_dir, anim_state) * sprite_size, ivec2(sprite_size))).center().alpha(alpha);
    }
};

struct World
{
    Worm worm;
    std::vector<Bug> bugs;
    Map map;
    ParticleController par;
    ivec2 camera_pos;

    float fade_in = 1;
    float fade_out = 0;

    std::optional<ivec2> hovered_tile_pos;
    bool hovered_tile_valid = false;
    bool any_tool_enabled = false;

    inline static Gui::StickyButton button_restart = {0, "[R] Restart level"};
    inline static Gui::Checkbox checkbox_pause = {2, 1, true, "[Space] Play/pause"};
    inline static Gui::Checkbox checkbox_halfspeed = {3, 4, false, "[Tab] Change game speed"};
    inline static Gui::ButtonList buttons_game_control = Gui::ButtonList(-screen_size.x / 2, -1, {&button_restart, &checkbox_halfspeed, &checkbox_pause});

    inline static Gui::RadioButtonList radiobuttons_tools;
    inline static Gui::RadioButton button_add_dirt = {6, radiobuttons_tools, "Create dirt"};
    inline static Gui::RadioButton button_add_bait = {8, radiobuttons_tools, "Create bait"};
    inline static Gui::RadioButton button_bomb = {5, radiobuttons_tools, "Destroy object"};
    inline static Gui::Button button_reverse = {7, "Reverse"};
    inline static Gui::ButtonList buttons_tools = Gui::ButtonList(screen_size.x / 2, 1, {&button_add_dirt, &button_add_bait, &button_reverse, &button_bomb});

    inline static const std::vector<Gui::ButtonList *> buttonlists = {&buttons_game_control, &buttons_tools};

    bool ShouldFadeOut() const
    {
        if (button_restart.was_clicked)
            return true; // Manual restart.
        if (worm.out_of_bounds)
            return true; // Out of bounds.
        if (worm.dead && worm.death_timer > 60)
            return true; // Death.
        if (!worm.segments.empty() && worm.segments.front() == map.data.exit_pos)
            return true; // Win condition.
        return false;
    }

    std::optional<int> ShouldChangeLevel()
    {
        if (ShouldFadeOut() && fade_out >= 1)
            return worm.dead || worm.out_of_bounds || button_restart.was_clicked ? map.level_index : map.level_index + 1;
        return {};
    }

    void Tick()
    {
        // The logic affected by game speed.
        for (int i = 0; i < (checkbox_halfspeed.enabled ? 1 : 2); i++)
        {
            // Worm.
            worm.Tick(map, par, !checkbox_pause.enabled);

            // Bugs.
            if (!checkbox_pause.enabled)
            {
                // Tick.
                for (Bug &bug : bugs)
                    bug.Tick(map, par);

                // Remove dead bugs.
                std::erase_if(bugs, [](const Bug &bug){return bug.death_timer >= bug.max_death_timer;});
            }

            // Worm-bug interaction.
            if (!checkbox_pause.enabled && !worm.dead)
            {
                for (ivec2 seg : worm.segments)
                {
                    if (std::any_of(bugs.begin(), bugs.end(), [&](const Bug &bug){return !bug.dead && bug.HitboxTilePos() == seg;}))
                    {
                        worm.dead = true;
                        worm.need_death_anim = true;
                        break;
                    }
                }
            }

            // Camera pos.
            camera_pos = map.tiles.size() * tile_size / 2;
            camera_pos.y += atlas.panel.size.y / 2;

            // Global timer.
            map.global_time++;
        }

        { // Listener pos.
            int separation = screen_size.x * 2;
            Audio::ListenerPosition(camera_pos.to_vec3(-separation));
            Audio::ListenerOrientation(fvec3(0,0,1), fvec3(0,-1,0));
            Audio::Source::DefaultRefDistance(separation);
        }

        // Particles.
        par.Tick();

        { // Gui.
            // Update button subscripts.
            button_add_dirt.subscript = map.data.num_dirt;
            button_add_bait.subscript = map.data.num_bait;
            button_bomb.subscript = map.data.num_bombs;
            button_reverse.subscript = map.data.num_reverses;

            { // Button blocking conditions.
                { // "Reverse worm" button.
                    bool &block = button_reverse.force_blocked;
                    block = false;

                    if (!block && (worm.dead || worm.out_of_bounds))
                        block = true; // Refuse to reverse if dead of out of bounds.
                    if (!block && std::any_of(worm.segments.begin(), worm.segments.end(), [&](ivec2 seg){return seg == map.data.worm_starting_pos;}))
                        block = true; // Refuse to reverse in the starting position.
                    if (!block && std::any_of(worm.segments.begin(), worm.segments.end(), [&](ivec2 seg){return seg == map.data.exit_pos;}))
                        block = true; // Refuse to reverse if touching the exit.
                }
            }

            // Button ticks.
            if (!ShouldFadeOut())
            {
                for (Gui::ButtonList *list : buttonlists)
                    list->Tick();
            }

            { // Extra keyboard and mouse controls.
                // Restart from keyboard.
                if (Input::Button(Input::r).pressed())
                    button_restart.was_clicked = true;
                // Pause from keyboard.
                if (Input::Button(Input::space).pressed())
                    checkbox_pause.enabled = !checkbox_pause.enabled;
                // Toggle speed from keyboard.
                if (Input::Button(Input::tab).pressed())
                    checkbox_halfspeed.enabled = !checkbox_halfspeed.enabled;

                // Unselect tool.
                if (mouse.right.pressed())
                {
                    for (size_t i = 0; i < radiobuttons_tools.con.Size(); i++)
                    {
                        if (radiobuttons_tools.con[i]->active)
                        {
                            radiobuttons_tools.con[i]->active = false;
                            Sounds::click(mouse.pos());
                        }
                    }
                }
            }
        }

        { // Using tools.
            // Check if any tool is enabled.
            any_tool_enabled = false;
            for (size_t i = 0; i < radiobuttons_tools.con.Size(); i++)
            {
                if (radiobuttons_tools.con[i]->active)
                {
                    any_tool_enabled = true;
                    break;
                }
            }

            // Detect hovered tile pos.
            hovered_tile_pos = {};
            hovered_tile_valid = false;
            if (any_tool_enabled && window.HasMouseFocus() && (mouse.pos().abs() < screen_size/2).all() && mouse.pos().y <= screen_size.y / 2 - atlas.panel.size.y)
            {
                hovered_tile_pos = div_ex(mouse.pos() + camera_pos, tile_size);
                if ((*hovered_tile_pos - 1 < 0).any() || (*hovered_tile_pos + 1 >= map.tiles.size()).any())
                    hovered_tile_pos = {}; // Refuse to interact with a tile on the map border, or outside.
            }

            // Specific tile-related tools.
            if (hovered_tile_pos)
            {
                if (button_bomb.active)
                {
                    // Bomb.
                    hovered_tile_valid = map.InfoAt(*hovered_tile_pos).bombable;

                    // Checking for a worm.
                    bool worm_hovered = false;
                    decltype(worm.segments)::iterator hovered_worm_segment;
                    if (!hovered_tile_valid && !worm.dead && (hovered_worm_segment = std::find(worm.segments.begin(), worm.segments.end(), *hovered_tile_pos)) != worm.segments.end())
                    {
                        hovered_tile_valid = true;
                        worm_hovered = true;
                    }

                    // Checking for a bug.
                    // Note that we don't remember the exact bug, which allows us to blow up more than one at a time.
                    bool bug_hovered = false;
                    if (!hovered_tile_valid && std::any_of(bugs.begin(), bugs.end(), [&](const Bug &bug){return bug.HitboxTilePos() == *hovered_tile_pos;}))
                    {
                        hovered_tile_valid = true;
                        bug_hovered = true;
                    }

                    if (hovered_tile_valid && mouse.left.pressed())
                    {
                        ivec2 pixel_pos = *hovered_tile_pos * tile_size + tile_size / 2;

                        if (worm_hovered)
                        {
                            // Shorten the worm.

                            if (worm.segments.back() == *hovered_tile_pos || (worm.segments.size() >= 2 && worm.segments[worm.segments.size() - 2] == *hovered_tile_pos))
                            {
                                // If one of the first two segments explodes, destroy the whole worm.
                                worm.dead = true;
                                worm.need_death_anim = true;
                                worm.no_sound_on_death_animation = true;
                            }
                            else
                            {
                                // Otherwise shrink the worm.
                                for (auto it = worm.segments.begin(); it != hovered_worm_segment + 1; it++)
                                    par.EffectBloodExplosion(*it * tile_size + tile_size / 2, fvec2(tile_size));
                                worm.segments.erase(worm.segments.begin(), hovered_worm_segment + 1);
                            }
                        }
                        else if (bug_hovered)
                        {
                            // Destroy bug.
                            for (auto it = bugs.begin(); it != bugs.end();)
                            {
                                bool is_hit = it->HitboxTilePos() == *hovered_tile_pos;
                                if (is_hit)
                                    it = bugs.erase(it);
                                else
                                    it++;
                            }
                            Sounds::bug_dies(pixel_pos);
                        }
                        else
                        {
                            // Destroy tile.
                            map.At(*hovered_tile_pos).type = Map::Tile::air;
                        }


                        par.EffectExplosion(pixel_pos, fvec2(tile_size));
                        Sounds::explosion(pixel_pos);
                        map.data.num_bombs--;
                    }
                }
                else if (button_add_dirt.active)
                {
                    // Add dirt.
                    hovered_tile_valid = map.At(*hovered_tile_pos).type == Map::Tile::air;
                    if (hovered_tile_valid && mouse.left.pressed())
                    {
                        map.At(*hovered_tile_pos).type = Map::Tile::dirt;
                        ivec2 pixel_pos = *hovered_tile_pos * tile_size + tile_size / 2;
                        par.EffectDirtMinor(pixel_pos, fvec2(tile_size), 6);
                        Sounds::dirt_placed(pixel_pos);

                        map.data.num_dirt--;
                    }
                }
                else if (button_add_bait.active)
                {
                    // Add bait.
                    hovered_tile_valid = map.At(*hovered_tile_pos).type == Map::Tile::air;
                    if (hovered_tile_valid && mouse.left.pressed())
                    {
                        map.At(*hovered_tile_pos).type = Map::Tile::bait;
                        ivec2 pixel_pos = *hovered_tile_pos * tile_size + tile_size / 2;
                        par.EffectBaitMinor(pixel_pos, fvec2(tile_size));
                        Sounds::bait(pixel_pos);

                        map.data.num_bait--;
                    }
                }
            }

            // Reverse worm.
            if (button_reverse.IsClicked())
            {
                map.data.num_reverses--;
                ivec2 midpoint = worm.segments[worm.segments.size() / 2];

                std::reverse(worm.segments.begin(), worm.segments.end());

                Sounds::reverse(midpoint * tile_size + tile_size / 2);

                for (size_t i = 0; i < radiobuttons_tools.con.Size(); i++)
                    radiobuttons_tools.con[i]->active = false;
            }
        }

        { // Fade.
            clamp_var_min(fade_in -= 0.055f);

            if (ShouldFadeOut())
                clamp_var_max(fade_out += 0.055);
        }
    }

    void Render() const
    {
        // Worm.
        worm.Draw(camera_pos);

        // Map.
        map.render(camera_pos, screen_size);

        // Bugs.
        for (const Bug &bug : bugs)
            bug.Render(camera_pos);

        // Particles.
        par.Render(camera_pos, screen_size);

        { // World gui.
            // Map border.
            if (any_tool_enabled)
            {
                constexpr int margin = 6;
                constexpr fvec3 color(0.8, 1, 0.5);
                constexpr float alpha = 0.15;
                constexpr float beta = 1;
                ivec2 border_pos(tile_size - margin);
                ivec2 border_size = (map.tiles.size() - 2) * tile_size + margin * 2;
                r.iquad(-camera_pos + border_pos, border_size with(y = 1)).color(color).alpha(alpha).beta(beta);
                r.iquad(-camera_pos + border_pos with(y++), border_size with(x = 1, y -= 2)).color(color).alpha(alpha).beta(beta);
                r.iquad(-camera_pos + border_pos with(y += border_size.y - 1), border_size with(y = 1)).color(color).alpha(alpha).beta(beta);
                r.iquad(-camera_pos + border_pos with(y++, x += border_size.x - 1), border_size with(x = 1, y -= 2)).color(color).alpha(alpha).beta(beta);
            }

            // Tile cursor.
            if (hovered_tile_pos)
                r.iquad(*hovered_tile_pos * tile_size + tile_size/2 - camera_pos, atlas.cursor.region(ivec2(atlas.cursor.size.y * !hovered_tile_valid, 0), ivec2(atlas.cursor.size.y))).alpha(0.9 + 0.1 * sin(window.Ticks() % 120 / 120.f * 2 * f_pi)).center();
        }

        // Fade.
        float fade = max(fade_in, fade_out);
        if (fade > 0)
            r.iquad(ivec2(0), screen_size).center().color(fvec3(0)).alpha(smoothstep(clamp(fade)));

        // Vignette.
        r.iquad(ivec2(0), atlas.vignette).center();

        { // Gui panel.
            // Panel background.
            r.iquad(ivec2(0, screen_size.y / 2), ivec2(screen_size.x, atlas.panel.size.y)).tex(atlas.panel.pos, atlas.panel.size).center(fvec2(0.5, atlas.panel.size.y));
            // Buttons.
            for (const Gui::ButtonList *list : buttonlists)
                list->Draw();

            // Level name
            if (float alpha = 1 - fade; alpha > 0)
            {
                constexpr fvec3 color = fvec3(131,185,195)/255;
                constexpr int half_spacing = 6;
                std::string index = FMT("#{}", map.level_index + 1);
                r.itext(ivec2(0, screen_size.y / 2 - atlas.panel.size.y / 2 - half_spacing * !map.data.level_name.empty()), Graphics::Text(Fonts::main, index)).color(color).alpha(alpha);
                if (!map.data.level_name.empty())
                    r.itext(ivec2(0, screen_size.y / 2 - atlas.panel.size.y / 2 + half_spacing), Graphics::Text(Fonts::main, FMT("\"{}\"", map.data.level_name))).color(color).alpha(alpha);
            }
        }

        // Mouse pointer.
        Draw::MousePointer();
    }
};

namespace States
{
    STRUCT( Game EXTENDS GameUtils::State::BasicState )
    {
        UNNAMED_MEMBERS()

        World w;
        Map map_backup;

        void LoadMap(int index)
        {
            if (map_backup.level_index != index)
            {
                map_backup = Map(index);

                { // Save progress.
                    Progress::Data progress = Progress::Load();
                    clamp_var_min(progress.level, index + 1);
                    Progress::Save(progress);
                }
            }
            ReloadCurrentMap();
        }

        void ReloadCurrentMap()
        {
            w.checkbox_pause.enabled = true;
            w.button_restart.was_clicked = false;
            for (size_t i = 0; i < w.radiobuttons_tools.con.Size(); i++)
            {
                w.radiobuttons_tools.con[i]->active = false;
                w.radiobuttons_tools.con[i]->subscript = 0;
            }

            w = World{};
            w.map = map_backup;

            for (int i = 0; i < w.map.data.worm_starting_len; i++)
                w.worm.segments.push_back(w.map.data.worm_starting_pos);

            for (const MapPoints::BugSpawn &spawn : w.map.data.bug_spawns)
                w.bugs.push_back(adjust(Bug{}, tile_pos = spawn.pos, dir = spawn.dir));

            // Tick, for a good measure.
            w.Tick();
        }

        Game()
        {
            Draw::InitAtlas();
            Music::StartOnce();
        }

        void Init(const std::string &params) override
        {
            int map_index = 0;
            if (!params.empty())
                map_index = Refl::FromString<int>(params);
            LoadMap(map_index);
        }

        void Tick(const GameUtils::State::NextStateSelector &next_state) override
        {
            (void)next_state;

            w.Tick();

            if (auto next_level = w.ShouldChangeLevel())
            {
                if (*next_level < Map::GetLevelCount())
                    LoadMap(*next_level);
                else
                    next_state.Set("Final");
            }

            // Toggle music.
            Music::ToggleIfHotkeyPressed();

            DebugInfo::Tick();
        }

        void Render() const override
        {
            Graphics::SetClearColor(fvec3(0));
            Graphics::Clear();

            r.BindShader();

            w.Render();

            DebugInfo::Render();
            r.Finish();
        }
    };

    STRUCT( Menu EXTENDS GameUtils::State::BasicState )
    {
        UNNAMED_MEMBERS()

        struct Button
        {
            Graphics::Text text;
            ivec2 pos{};
            ivec2 size{};
            fvec2 align{};

            Button() {}

            Button(ivec2 pos, fvec2 align, std::string_view str)
                : pos(pos), align(align)
            {
                SetText(str);
            }

            void SetText(std::string_view str)
            {
                text = Graphics::Text(Fonts::main, str);
                size = text.ComputeStats().size;
            }

            bool IsHovered() const
            {
                constexpr ivec2 padding(4, 4);

                fvec2 offset = -size * (align * 0.5 + 0.5);
                return (mouse.pos() >= pos + offset - padding).all() && (mouse.pos() < pos + offset + size + padding).all();
            }

            bool IsClicked() const
            {
                return mouse.left.released() && IsHovered();
            }

            void Draw() const
            {
                fvec3 color = IsHovered() ? (mouse.left.down() ? fvec3(109, 212, 12) / 255 : fvec3(247, 251, 14) / 255) : fvec3(0.95);
                r.itext(pos, text).align(align).color(color);
            }
        };

        struct Continue
        {
            Button b_continue;
            Button b_minus, b_plus;
        };
        std::optional<Continue> cont;
        Button b_new_game, b_quit;

        Progress::Data progress;

        int cont_level = -1; // 0-based.

        float fade = 0;
        std::optional<int> target_level_index; // 0-based.

        static std::string GetContinueString(int level_index)
        {
            return FMT("Continue from level #{}", level_index);
        }

        Menu()
        {
            Draw::InitAtlas();
            progress = Progress::Load();
            clamp_var_max(progress.level, Map::GetLevelCount());
            cont_level = progress.level;

            constexpr int button_spacing = 24, plus_minus_spacing = 16;

            // Note that we use a hardcoded level index here, because we don't want the real one to affect positions.
            int continue_string_w = Graphics::Text(Fonts::main, GetContinueString(42/*sic*/)).ComputeStats().size.x;

            // The X pos of the +/- buttons.
            int plus_minus_pos_x = continue_string_w / 2 + plus_minus_spacing;

            int pos_y = 0;
            if (cont_level > 0)
            {
                // "Continue" buttons, if needed.
                cont.emplace();
                cont->b_continue = Button(ivec2(-continue_string_w/2, pos_y), fvec2(-1, 0), GetContinueString(cont_level));
                cont->b_minus = Button(ivec2(plus_minus_pos_x, pos_y), fvec2(0), "<");
                cont->b_plus = Button(ivec2(plus_minus_pos_x + plus_minus_spacing, pos_y), fvec2(0), ">");
            }
            pos_y += button_spacing;
            b_new_game = Button(ivec2(0, pos_y), fvec2(0), "New game");
            pos_y += button_spacing;
            b_quit = Button(ivec2(0, pos_y), fvec2(0), "Quit");
        }

        void Tick(const GameUtils::State::NextStateSelector &next_state) override
        {
            (void)next_state;

            if (!target_level_index)
            {
                if (cont)
                {
                    if (cont->b_continue.IsClicked())
                    {
                        target_level_index = cont_level - 1;
                    }
                    if (cont->b_minus.IsClicked())
                    {
                        clamp_var_min(--cont_level, 1);
                        cont->b_continue.SetText(GetContinueString(cont_level));
                    }
                    if (cont->b_plus.IsClicked())
                    {
                        clamp_var_max(++cont_level, progress.level);
                        cont->b_continue.SetText(GetContinueString(cont_level));
                    }
                }

                if (b_new_game.IsClicked())
                    target_level_index = 0;

                if (b_quit.IsClicked())
                    Program::Exit();
            }
            else
            {
                if (fade >= 1)
                    next_state.Set("Game", Refl::ToString(target_level_index.value_or(0)));
                clamp_var_max(fade += 0.025);
            }

            // Toggle music.
            Music::ToggleIfHotkeyPressed();

            DebugInfo::Tick();
        }

        void Render() const override
        {
            Graphics::SetClearColor(fvec3(0));
            Graphics::Clear();

            r.BindShader();

            if (cont)
            {
                cont->b_continue.Draw();
                cont->b_plus.Draw();
                cont->b_minus.Draw();
            }
            b_new_game.Draw();
            b_quit.Draw();

            constexpr std::string_view info_string = "Press [M] at any time to mute music.\nA game by HolyBlackCat for Ludum dare #48, Apr 2021.";
            r.itext(ivec2(0, screen_size.y / 2), Graphics::Text(Fonts::main, info_string)).color(fvec3(1)).alpha(0.3).align(ivec2(0,1));

            Draw::MousePointer();

            if (fade > 0)
                r.iquad(ivec2(0), screen_size).center().color(fvec3(0)).alpha(smoothstep(fade));

            DebugInfo::Render();
            r.Finish();
        }
    };

    STRUCT( Final EXTENDS GameUtils::State::BasicState )
    {
        UNNAMED_MEMBERS()

        int time = 0;

        int worm_segments = 8;

        Final()
        {
            Draw::InitAtlas();
            Music::StartOnce();
        }

        void Tick(const GameUtils::State::NextStateSelector &next_state) override
        {
            (void)next_state;
            if (time > 800)
                Program::Exit();

            time++;

            // Toggle music.
            Music::ToggleIfHotkeyPressed();

            DebugInfo::Tick();
        }

        void Render() const override
        {
            Graphics::SetClearColor(fvec3(0));
            Graphics::Clear();

            r.BindShader();

            constexpr int y_offset = 12;

            { // Worm.
                ivec2 pos(-screen_size.x / 2 - 60 + time, y_offset);

                // Head.
                r.iquad(pos, atlas.worm.region(ivec2(1,1) * tile_size, ivec2(tile_size))).center();
                pos.x -= tile_size;

                // Body.
                for (int i = 0; i < worm_segments; i++)
                {
                    r.iquad(pos, atlas.worm.region(ivec2(0,0) * tile_size, ivec2(tile_size))).center();
                    pos.x -= tile_size;
                }

                // Tail.
                r.iquad(pos, atlas.worm.region(ivec2(2,2) * tile_size, ivec2(tile_size))).center();
            }

            constexpr int t1 = 200, t2 = 320;
            constexpr int text_anim_len = 60;
            r.itext(ivec2(0,y_offset - 32), Graphics::Text(Fonts::main, "That's all!")).alpha(smoothstep(clamp(Math::linear_mapping<float>(t1, t1 + text_anim_len, 0, 1)(time))));
            r.itext(ivec2(0,y_offset + 32), Graphics::Text(Fonts::main, "Thanks for playing my game")).alpha(smoothstep(clamp(Math::linear_mapping<float>(t2, t2 + text_anim_len, 0, 1)(time))));

            DebugInfo::Render();
            r.Finish();
        }
    };
}
