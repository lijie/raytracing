// 啊... GLSL 写不了 interface
// 啊... 也没有 void* 或者 any 这样的类型
// 啊... 数组作为函数参数时必须定长
// :(

// 射线
struct ray {
    vec3 origin;    // 原点
    vec3 direction; // 方向, 单位向量
};

ray ray_new(vec3 origin, vec3 direction) {
    ray r;
    r.origin = origin;
    r.direction = direction;
    return r;
}

vec3 ray_point_at(ray r, float t) {
    return r.origin + r.direction * t;
}

struct camera {
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 origin;
};

camera camera_new() {
    camera c;
    c.lower_left_corner = vec3(-2.0, -1.0, -1.0);
    c.horizontal = vec3(4.0, 0.0, 0.0);
    c.vertical = vec3(0.0, 2.0, 0.0);
    c.origin = vec3(0.0, 0.0, 0.0);
    return c;
}

ray camera_get_ray(camera c, vec2 uv) {
    vec3 direction = normalize(c.lower_left_corner + uv.x * c.horizontal + uv.y * c.vertical - c.origin);
    return ray_new(c.origin, direction);
}

struct hit_record {
    vec3 p;         // 交点
    vec3 normal;    // 交点处的法线
    float t;        // 交点与射线起点的距离
    float front;    // 交点是在几何表面, 还是在里面
};

void hit_record_set_face_normal(out hit_record rec, ray r, vec3 out_normal) {
    // 如果在表面
    // 则射线(来自观察者)与法线的点积应该<0
    rec.front = -dot(r.direction, out_normal);
    if (rec.front > 0.0)
        rec.normal = out_normal;
    else
        rec.normal = -out_normal;
}

struct sphere {
    vec3 center;
    float radius;
};

sphere sphere_new(vec3 center, float radius) {
    sphere sp;
    sp.center = center;
    sp.radius = radius;
    return sp;
}

float sphere_hit(sphere sp, ray r, float t_min, float t_max, out hit_record rec) {
    // 计算一元二次方程
    // b^2 - 4ac
    float a = dot(r.direction, r.direction);
    float b = 2.0 * dot(r.origin - sp.center, r.direction);
    float c = dot((r.origin - sp.center), (r.origin - sp.center)) - sp.radius * sp.radius;
    float discriminant = b * b - 4.0 * a * c;
    if (discriminant < 0.0)
        return -1.0;
    
    float root = (-(b) - sqrt(discriminant)) / (2.0 * a);

    if (root <= t_min || root >= t_max) {
        root = (-(b) + sqrt(discriminant)) / (2.0 * a);
        if (root <= t_min || root >= t_max) {
            // 无解
            return 0.0;
        }
    }

    // 有效解
    rec.t = root;
    rec.p = ray_point_at(r, rec.t);
    vec3 out_normal = normalize(rec.p - sp.center);
    hit_record_set_face_normal(rec, r, out_normal);
    return 1.0;
}

sphere[100] world_sphere_array;
int world_sphere_count = 0;

float world_hit(ray r, float t_min, float t_max, out hit_record rec) {
    float closest = t_max;
    float hit = 0.0;
    for (int i = 0; i < world_sphere_count; i++) {
        sphere sp = world_sphere_array[i];
        hit = sphere_hit(sp, r, t_min, closest, rec);
        if (hit > 0.0) {
            closest = rec.t;
        }
    }
    return t_max - closest;
}

void world_add(sphere sp) {
    world_sphere_array[world_sphere_count] = sp;
    world_sphere_count++;
}

// 自顶向下返回一个蓝白混合的颜色
vec3 ray_color(ray r) {
    // sphere sp = sphere_new(vec3(0, 0, -1), 0.5);
    hit_record rec;
    // float t = sphere_hit(sp, r, 0.0, 10000.0, rec);
    float t = world_hit(r, 0.0, 100000.0, rec);
    if (t > 0.0) {
        // scale to [0, 1]
        return 0.5 * (rec.normal + 1.0);
    }
    // scale t to [0, 1]
    t = 0.5 * (r.direction.y + 1.0);
    return mix(vec3(1.0, 1.0, 1.0), vec3(0.5, 0.7, 1.0), t);
}

void main()
{
    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    camera c = camera_new();
    ray r = camera_get_ray(c, uv);
    world_add(sphere_new(vec3(0.0, 0.0,-1.0), 0.5));
    world_add(sphere_new(vec3(0.0, -100.5,-1.0), 100.0));
    vec3 color = ray_color(r);
    gl_FragColor = vec4(color, 1);
}
