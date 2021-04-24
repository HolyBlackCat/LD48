#include "game/main.h"

SIMPLE_STRUCT( Atlas
    DECL(Graphics::TextureAtlas::Region)
        worm
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
        for (std::size_t i = 0; i < segments.size(); i++)
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

        Game()
        {
            static bool once = true;
            if (once)
            {
                once = false;
                texture_atlas.InitRegions(atlas, ".png");
            }

            Refl::FromString(worm.segments, "[(0,0),(1,0),(2,0),(3,0),(4,0),(5,0),(5,-1),(5,-2),(5,-3),(5,-4),(4,-4),(3,-4),(2,-4),(1,-4),(0,-4),(-1,-4),(-2,-4),(-3,-4),(-4,-4),(-5,-4),(-5,-3),(-5,-2),(-5,-1),(-5,0),(-5,1),(-5,2),(-5,3),(-4,3),(-3,3),(-2,3),(-1,3),(0,3),(1,3),(2,3),(3,3),(4,3),(5,3),(6,3),(7,3),(8,3),(9,3),(9,2),(9,1),(9,0),(9,-1),(9,-2),(9,-3),(9,-4),(10,-4),(11,-4),(12,-4),(13,-4),(14,-4),(15,-4),(16,-4),(17,-4),(18,-4),(19,-4),(20,-4),(20,-3),(20,-2),(20,-1),(20,0),(20,1),(20,2),(20,3),(20,4),(19,4),(18,4),(17,4),(16,4),(15,4),(14,4),(13,4),(12,4),(12,3),(12,2),(12,1),(12,0),(12,-1),(13,-1),(14,-1),(15,-1),(16,-1),(17,-1)]");
        }

        void Tick(const GameUtils::State::NextStateSelector &next_state) override
        {
            (void)next_state;

            ivec2 dir(b_right.repeated() - b_left.repeated(), b_down.repeated() - b_up.repeated());
            if (abs(dir.x) != abs(dir.y))
                worm.segments.push_back(worm.segments.empty() ? ivec2(0,0) : worm.segments.back() + dir);

            if (ImGui::Button("Save"))
            {
                SDL_SetClipboardText(Refl::ToString(worm.segments).c_str());
                std::cout << Refl::ToString(worm.segments) << '\n';
            }

            worm.crawl_offset += Input::Button(Input::mouse_wheel_up).pressed() - Input::Button(Input::mouse_wheel_down).pressed();

            if (Input::Button(Input::space).down())
                worm.Crawl(1);
        }

        void Render() const override
        {
            Graphics::SetClearColor(fvec3(0));
            Graphics::Clear();

            r.BindShader();

            worm.Draw(mouse.pos());

            r.itext(-screen_size / 2, Graphics::Text(Fonts::main, STR((worm.crawl_offset)))).align(ivec2(-1));

            r.Finish();
        }
    };
}
