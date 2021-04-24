#include "game/main.h"
#include "gameutils/tiled_map.h"
#include "utils/json.h"

SIMPLE_STRUCT( Atlas
    DECL(Graphics::TextureAtlas::Region)
        worm, tiles
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
}

struct Worm
{
    std::vector<ivec2> segments; // This is sorted tail-to-head.
    int crawl_offset = 0;

    void Draw(ivec2 offset) const
    {
        std::optional<ivec2> prev_delta, next_delta;
        for (size_t i = 0; i < segments.size(); i++)
        {
            prev_delta = next_delta;
            next_delta = i == segments.size() - 1 ? std::nullopt : std::optional{segments[i+1] - segments[i]};

            Graphics::TextureAtlas::Region region = atlas.worm.region(ivec2(0), ivec2(tile_size));

            ivec2 base_pos = offset + segments[i] * tile_size;

            if (!prev_delta && !next_delta)
            {
                // This shouldn't happen, draw a placeholder.
                r.iquad(base_pos, region);
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

                Draw::ShiftedTile(base_pos, {region, nullptr, &straight_reg, (*next_delta < 0).any(), true}, *next_delta * crawl_offset);
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

                Draw::ShiftedTile(base_pos, {region, &straight_reg, nullptr, (*prev_delta < 0).any(), true}, *prev_delta * crawl_offset);
            }
            else if (*prev_delta == *next_delta)
            {
                // Straight piece.
                if (next_delta->x != 0)
                    region = region.region(ivec2(0,0) * tile_size, ivec2(tile_size));
                else
                    region = region.region(ivec2(1,0) * tile_size, ivec2(tile_size));

                Draw::ShiftedTile(base_pos, {region, i == 1 ? nullptr : &region, i == segments.size() - 2 ? nullptr : &region, (*prev_delta < 0).any() }, *next_delta * crawl_offset);
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
                r.iquad(base_pos, region);
            }
        }
        std::cout << '\n';
    }

    ivec2 ChooseDirection() const
    {
        return ivec2(Input::Button(Input::d).down() - Input::Button(Input::a).down(), Input::Button(Input::s).down() - Input::Button(Input::w).down());
    }

    void Crawl(int delta)
    {
        if (segments.size() < 2)
            return;


        // We crawl pixel by pixel, just in case.
        while (delta-- > 0)
        {
            crawl_offset++;
            if (crawl_offset >= tile_size / 2)
            {
                crawl_offset -= tile_size;
                segments.front() = segments.back() + ChooseDirection();
                std::rotate(segments.begin(), segments.begin() + 1, segments.end());
            }
        }
    }
};

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

        r.fquad(pos, fvec2(params.size)).center().color(params.color).alpha(params.alpha).beta(params.beta);
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

    void EffectDirtDestroyed(fvec2 center_pos, fvec2 area, int n = 15)
    {
        while (n-- > 0)
        {
            float p = 0 <= randomize.real() <= 1;

            particles.push_back(adjust(Particle{}
                , pos = center_pos + fvec2(-area.x/2 <= randomize.real() <= area.x/2, -area.y/2 <= randomize.real() <= area.y/2)
                , vel = fvec2::dir(randomize.angle(), 0 <= randomize.real() <= 1)
                , acc = fvec2(0,0.05)
                , max_time = 30 <= randomize.integer() <= 60
                , damp = mix(p, 0.01, 0.03)
                , p.size = mix(p, 8, 2)
                , p.color = (std::array{fvec3(139, 90, 60) / 255, fvec3(101, 62, 41) / 255, fvec3(159, 111, 56) / 255}[randomize.integer() < 3] * float(0.9 <= randomize.real() <= 1.1))
                // , p2 = _object_.p
                // , p2->t = _object_.max_time / 2
                , p1 = _object_.p
                , p1->t = _object_.max_time
                , p1->size = 0
            ));
        }
    }
};

class Map
{
  public:
    // Tile types.
    enum class Tile
    {
        air,
        dirt,
        pipe,
        _count,
    };

    enum class TileRenderer {simple, fancy, pipe};

    struct TileInfo
    {
        bool solid = true;
        TileRenderer renderer{};
        int tile_index = -1;
    };

    // Tile properties table.
    inline static const TileInfo tile_info[] =
    {
        /* air  */ {.solid = false, .renderer = TileRenderer::simple, .tile_index = -1},
        /* dirt */ {.solid =  true, .renderer = TileRenderer::fancy , .tile_index =  0},
        /* pipe */ {.solid =  true, .renderer = TileRenderer::pipe  , .tile_index =  1},
    };
    static_assert(std::size(tile_info) == size_t(Tile::_count));

    // Tile merging rules.
    static bool ShouldMergeWith(Tile a, Tile b)
    {
        if (a == b)
            return true;
        return false;
    }


    struct TileData
    {
        Tile type = Tile::air;
    };

    Array2D<TileData> tiles;
    int level_index = -1;

    Array2D<uint8_t> noise;

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
        }
        catch (std::exception &e)
        {
            Program::Error(STR("While loading map #", (new_level_index+1), ":\n", (e.what())));
        }
    }

    const TileData &At(ivec2 pos) const
    {
        return tiles.clamped_at(pos);
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

                    auto corner = [&](ivec2 sub)
                    {
                        ivec2 tile_offset = sub * 2 - 1;
                        ivec2 pixel_offset = sub * (tile_size / 2);

                        bool merge_h = merge_full || ShouldMergeWith(type, At(tile_pos + tile_offset with(y = 0)).type);
                        bool merge_v = merge_full || ShouldMergeWith(type, At(tile_pos + tile_offset with(x = 0)).type);
                        bool merge_d = merge_full || ShouldMergeWith(type, At(tile_pos + tile_offset            ).type);

                        Graphics::TextureAtlas::Region base_region = atlas.tiles.region(ivec2(1, 2 * info.tile_index) * tile_size, ivec2(tile_size));
                        Graphics::TextureAtlas::Region region;

                        if (merge_h && merge_v && merge_d)
                        {
                            region = base_region;
                            if (merge_full)
                            {
                                int variants[] = {0,0,0,0,0,1,1,1,2,3};
                                int var = variants[Noise(tile_pos) % std::size(variants)];
                                region.pos += bvec2(var & 1, var & 2) * tile_size;
                            }
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

                    for (ivec2 sub : vector_range(ivec2(2)))
                        corner(sub);
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
                    r.iquad(pixel_pos, atlas.tiles.region((ivec2(1, 2 * info.tile_index) + variant) * tile_size, ivec2(tile_size)));
                }
                break;
            }
        }
    }
};

namespace States
{
    STRUCT( Game EXTENDS GameUtils::State::BasicState )
    {
        UNNAMED_MEMBERS()

        Worm worm;

        Input::Button b_up = Input::up;
        Input::Button b_down = Input::down;
        Input::Button b_left = Input::left;
        Input::Button b_right = Input::right;

        Map map;

        ParticleController par;

        Game()
        {
            static bool once = true;
            if (once)
            {
                once = false;
                texture_atlas.InitRegions(atlas, ".png");
            }

            Refl::FromString(worm.segments, "[(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(5,-1),(5,-2),(5,-3),(5,-4),(4,-4),(3,-4),(2,-4),(1,-4),(0,-4),(-1,-4),(-2,-4),(-3,-4),(-4,-4),(-5,-4),(-5,-3),(-5,-2),(-5,-1),(-5,0),(-5,1),(-5,2),(-5,3),(-4,3),(-3,3),(-2,3),(-1,3),(0,3),(1,3),(2,3),(3,3),(4,3),(5,3),(6,3),(7,3),(8,3),(9,3),(9,2),(9,1),(9,0),(9,-1),(9,-2),(9,-3),(9,-4),(10,-4),(11,-4),(12,-4),(13,-4),(14,-4),(15,-4),(16,-4),(17,-4),(18,-4),(19,-4),(20,-4),(20,-3),(20,-2),(20,-1),(20,0),(20,1),(20,2),(20,3),(20,4),(19,4),(18,4),(17,4),(16,4),(15,4),(14,4),(13,4),(12,4),(12,3),(12,2),(12,1),(12,0),(12,-1),(13,-1),(14,-1),(15,-1),(16,-1),(17,-1)]");

            map = Map(0);
        }

        void Tick(const GameUtils::State::NextStateSelector &next_state) override
        {
            (void)next_state;

            { // Worm
                ivec2 dir(b_right.repeated() - b_left.repeated(), b_down.repeated() - b_up.repeated());
                if (abs(dir.x) != abs(dir.y))
                    worm.segments.push_back(worm.segments.empty() ? ivec2(0,0) : worm.segments.back() + dir);

                worm.crawl_offset += Input::Button(Input::mouse_wheel_up).pressed() - Input::Button(Input::mouse_wheel_down).pressed();

                if (Input::Button(Input::space).down())
                    worm.Crawl(1);
            }

            { // Particles
                par.Tick();

                if (mouse.left.released())
                    par.EffectDirtDestroyed(mouse.pos(), fvec2(16));
            }
        }

        void Render() const override
        {
            Graphics::SetClearColor(fvec3(0));
            Graphics::Clear();

            r.BindShader();

            worm.Draw(mouse.pos());

            r.itext(-screen_size / 2, Graphics::Text(Fonts::main, STR((worm.crawl_offset)))).align(ivec2(-1));

            ivec2 camera_pos = mouse.pos();

            map.render(camera_pos, screen_size);
            par.Render(camera_pos, screen_size);

            r.Finish();
        }
    };
}
