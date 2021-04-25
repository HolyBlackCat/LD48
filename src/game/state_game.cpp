#include "game/main.h"
#include "gameutils/tiled_map.h"
#include "signals/connection.h"
#include "utils/json.h"

SIMPLE_STRUCT( Atlas
    DECL(Graphics::TextureAtlas::Region)
        worm, tiles, vignette, panel, button, button_icons
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
    ADD_SOUND(landing, 0.2)
    ADD_SOUND(click, 0.03)
    ADD_SOUND(player_dies, 0.1)

    #undef ADD_SOUND
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

        bool IsHovered() const
        {
            return (mouse.pos() >= pos).all() && (mouse.pos() < pos + Size()).all();
        }

        bool IsClicked() const
        {
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

        virtual void Draw() const
        {
            int state = IsHovered() ? (mouse.left.down() ? 2 : 1) : 0;
            r.iquad(pos, atlas.button.region(ivec2(SideDecoW() + Button::Size() * state, 0), ivec2(Button::Size())));
            r.iquad(pos, atlas.button_icons.region(ivec2(Button::Size() * GetIconIndex(), 0), ivec2(Button::Size()))).alpha(state == 2 ? 0.75 : 1);
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
        void Click() override
        {
            Button::Click();
            was_clicked = true;
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

    void EffectDirtDestroyed(fvec2 center_pos, fvec2 area, int n = 15)
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
                fvec2::dir(randomize.angle(), 0 <= randomize.real() <= 2),
                2, 16, 45 <= randomize.integer() <= 90)
            );
        }
    }
};

namespace MapPoints
{
    struct State
    {
        bool worm_placed = false;
        ivec2 worm_starting_pos{};
        ivec2 worm_starting_dir{};
        int worm_starting_len = 4;

        bool exit_placed = false;
        ivec2 exit_pos{};
        ivec2 exit_dir{};
    };

    STRUCT( Base POLYMORPHIC )
    {
        UNNAMED_MEMBERS()
        virtual void Process(State &state, ivec2 tile_pos) = 0;
    };

    STRUCT( Worm EXTENDS Base )
    {
        MEMBERS(
            DECL(int) len
            DECL(ivec2) dir
        )

        void Process(State &state, ivec2 tile_pos) override
        {
            if (state.worm_placed)
                Program::Error("More than one worm.");
            state.worm_starting_pos = tile_pos;
            state.worm_starting_dir = dir;
            state.worm_starting_len = len;
        }
    };

    STRUCT( Exit EXTENDS Base )
    {
        MEMBERS(
            DECL(ivec2) dir
        )

        void Process(State &state, ivec2 tile_pos) override
        {
            if (state.worm_placed)
                Program::Error("More than one exit.");
            state.exit_placed = true;
            state.exit_pos = tile_pos;
            state.exit_dir = dir;
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
    };

    // Tile properties table.
    inline static const TileInfo tile_info[] =
    {
        /* air   */ {.renderer = TileRenderer::simple, .tile_index = -1, .breakable = false, .path_priority = 99, .kills = false},
        /* dirt  */ {.renderer = TileRenderer::fancy , .tile_index =  0, .breakable = true , .path_priority =  9, .kills = false},
        /* pipe  */ {.renderer = TileRenderer::pipe  , .tile_index =  1, .breakable = false, .path_priority =  0, .kills = false},
        /* wall  */ {.renderer = TileRenderer::pipe  , .tile_index =  2, .breakable = false, .path_priority =  0, .kills = false},
        /* spike */ {.renderer = TileRenderer::pipe  , .tile_index =  3, .breakable = false, .path_priority =  9, .kills = true },
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

struct Worm
{
    std::vector<ivec2> segments; // This is sorted tail-to-head.
    int crawl_offset = 0;
    float fall_offset = 0;
    float fall_speed = 0;
    bool dead = false;
    int death_timer = 0;

    ivec2 ChooseDirection(const Map &map) const
    {
        if (segments.size() < 2)
            return ivec2();
        ivec2 head = segments.back();
        ivec2 pre_head = segments[segments.size() - 2];
        ivec2 current_dir = head - pre_head;
        if (!current_dir)
            current_dir = map.data.worm_starting_dir;

        auto GetPathPriority = [&](ivec2 tile_pos) -> int
        {
            int ret = 0;

            if (std::find(segments.begin(), segments.end(), tile_pos) != segments.end())
                ret = 0; // Refuse to touch your own body.
            else if (tile_pos == map.data.exit_pos && tile_pos - segments.back() == -map.data.exit_dir)
                ret = 98; // Go the exit if possible.
            else
                ret = map.InfoAt(tile_pos).path_priority;

            return ret;
        };

        auto p_fwd = GetPathPriority(segments.back() + current_dir);
        auto p_left = GetPathPriority(segments.back() + current_dir.rot90(-1));
        auto p_right = GetPathPriority(segments.back() + current_dir.rot90(1));

        // Stop if nowhere to go.
        if (max(p_fwd, p_left, p_right) <= 0)
            return ivec2();

        // If the current dir is at least as good as the alternatives, keep it.
        if (p_fwd >= p_left && p_fwd >= p_right)
            return current_dir;

        // Pick the best of the two directions.
        if (p_left > p_right)
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
                    break;
                }
            }

            // Death animation.
            if (dead)
            {
                for (ivec2 seg : segments)
                    par.EffectBloodExplosion(seg * tile_size + tile_size/2, fvec2(tile_size));

                Sounds::player_dies(segments[segments.size() / 2] * tile_size + tile_size / 2);
            }
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
                    bool at_exit = map.data.exit_pos == segments.back();
                    ivec2 chosen_dir = at_exit ? ivec2() : ChooseDirection(map);

                    if (chosen_dir || at_exit)
                    {
                        crawl_offset -= tile_size;
                        segments.front() = segments.back() + chosen_dir;
                        std::rotate(segments.begin(), segments.begin() + 1, segments.end());

                        { // Destroy tiles.
                            if (map.InfoAt(segments.back()).breakable)
                            {
                                map.At(segments.back()).type = Map::Tile::air;
                                ivec2 pixel_pos = segments.back() * tile_size + tile_size / 2;
                                par.EffectDirtDestroyed(pixel_pos, fvec2(tile_size));
                                Sounds::dirt_breaks(pixel_pos, 0.5);
                            }
                        }
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
        float alpha = clamp_min(1 - death_timer / 45.f);
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

                Draw::ShiftedTile(base_pos, {region, &straight_reg, nullptr, (*prev_delta < 0).any(), true}, *prev_delta * crawl_offset, alpha);
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

namespace Controls
{
    bool PausePressed()
    {
        return mouse.right.pressed() || Input::Button(Input::space).pressed();
    }
    bool ToggleSpeedPressed()
    {
        return Input::Button(Input::tab).pressed();
    }
}

struct World
{
    Worm worm;
    Map map;
    ParticleController par;
    ivec2 camera_pos;

    float fade_in = 1;
    float fade_out = 0;

    inline static Gui::StickyButton button_restart = {0, "[R] Restart level"};
    inline static Gui::Checkbox checkbox_pause = {2, 1, true, "[Space] Play/pause"};
    inline static Gui::Checkbox checkbox_halfspeed = {3, 4, false, "[Tab] Change game speed"};
    inline static Gui::ButtonList buttons_mid = Gui::ButtonList(-screen_size.x / 2, -1, {&button_restart, &checkbox_halfspeed, &checkbox_pause});

    bool ShouldFadeOut() const
    {
        if (button_restart.was_clicked)
            return true; // Manual restart.
        if (worm.dead && worm.death_timer > 60)
            return true; // Death.
        if (!worm.segments.empty() && worm.segments.front() == map.data.exit_pos)
            return true; // Win condition.
        return false;
    }

    std::optional<int> ShouldChangeLevel()
    {
        if (ShouldFadeOut() && fade_out >= 1)
            return worm.dead || button_restart.was_clicked ? map.level_index : map.level_index + 1;
        return {};
    }

    void Tick()
    {
        // The logic affected by game speed.
        for (int i = 0; i < (checkbox_halfspeed.enabled ? 1 : 2); i++)
        {
            { // Worm
                worm.Tick(map, par, !checkbox_pause.enabled);
            }

            // Camera pos.
            camera_pos = map.tiles.size() * tile_size / 2;

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
            if (!ShouldFadeOut())
            {
                buttons_mid.Tick();
            }

            if (Controls::PausePressed())
                checkbox_pause.Click();

            if (Controls::ToggleSpeedPressed())
                checkbox_halfspeed.Click();
        }

        { // Fade.
            clamp_var_min(fade_in -= 0.035f);

            if (ShouldFadeOut())
                clamp_var_max(fade_out += 0.035f);
        }
    }

    void Render() const
    {
        worm.Draw(camera_pos);
        map.render(camera_pos, screen_size);
        par.Render(camera_pos, screen_size);

        { // Gui.
            r.iquad(ivec2(0, screen_size.y / 2), ivec2(screen_size.x, atlas.panel.size.y)).tex(atlas.panel.pos, atlas.panel.size).center(fvec2(0.5, atlas.panel.size.y));
            buttons_mid.Draw();
        }

        // Fade.
        if (float fade = max(fade_in, fade_out); fade > 0)
            r.iquad(ivec2(0), screen_size).center().color(fvec3(0)).alpha(fade);

        // Vignette.
        r.iquad(ivec2(0), atlas.vignette).center();
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
                map_backup = Map(index);
            ReloadCurrentMap();
        }

        void ReloadCurrentMap()
        {
            w.checkbox_pause.enabled = true;
            w.button_restart.was_clicked = false;
            w = World{};
            w.map = map_backup;
            for (int i = 0; i < w.map.data.worm_starting_len; i++)
                w.worm.segments.push_back(w.map.data.worm_starting_pos);
        }

        Game()
        {
            static bool once = true;
            if (once)
            {
                once = false;
                texture_atlas.InitRegions(atlas, ".png");
            }

            LoadMap(0);
        }

        void Tick(const GameUtils::State::NextStateSelector &next_state) override
        {
            (void)next_state;

            w.Tick();

            if (auto next_level = w.ShouldChangeLevel())
                LoadMap(*next_level);
        }

        void Render() const override
        {
            Graphics::SetClearColor(fvec3(0));
            Graphics::Clear();

            r.BindShader();

            w.Render();

            r.Finish();
        }
    };
}
