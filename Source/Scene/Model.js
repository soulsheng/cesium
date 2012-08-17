/*global define*/
define([
        '../Core/DeveloperError',
        '../Renderer/BufferUsage',
        '../Core/ComponentDatatype',
        '../Core/IndexDatatype',
        '../Core/PrimitiveType',
        '../Renderer/CullFace',
        '../Core/Matrix4'
    ], function(
        DeveloperError,
        BufferUsage,
        ComponentDatatype,
        IndexDatatype,
        PrimitiveType,
        CullFace,
        Matrix4
        ) {
    "use strict";

    var Model = function(json) {
        this._data = json;

        this._sp = undefined;
        this._rs = undefined;

        this._vas = undefined;

        this.scale = 1.0;

        this.modelMatrix = Matrix4.IDENTITY;

        var that = this;
        this._drawUniforms = {
            u_scale : function() {
                return that.scale;
            },
            u_model : function() {
                return that.modelMatrix;
            }
        };
    };

    function getGeometry(context, geometry) {
        if (typeof geometry === 'undefined' || typeof geometry.positions === 'undefined' ||
                (typeof geometry.indices === 'undefined' && typeof geometry.nodes === 'undefined')) {
            return undefined;
        }

        var datatype = ComponentDatatype.FLOAT;
        var usage = BufferUsage.STATIC_DRAW;
        var positions = datatype.toTypedArray(geometry.positions);
        var pb = context.createVertexBuffer(positions, usage);

        var attributes = [{
            index : 0,
            vertexBuffer : pb,
            componentDatatype : datatype,
            componentsPerAttribute : 3
        }];

        if (typeof geometry.normals !== 'undefined') {
            var normals = datatype.toTypedArray(geometry.normals);
            var nb = context.createVertexBuffer(normals, usage);

            attributes.push({
                index : 1,
                vertexBuffer : nb,
                componentDatatype : datatype,
                componentsPerAttribute : 3
            });
        }

        var ib;
        var va;

        if (typeof geometry.indices !== 'undefined') {
            ib = context.createIndexBuffer(new Uint16Array(geometry.indices), usage, IndexDatatype.UNSIGNED_SHORT);
            va = context.createVertexArray(attributes, ib);

            return [{
                id : geometry.coreId,
                type : PrimitiveType.TRIANGLES,
                vertexArray : va
            }];
        } else if (typeof geometry.nodes !== 'undefined') {
            var vas = [];
            var nodes = geometry.nodes;
            var length = nodes.length;
            for (var i = 0; i < length; ++i) {
                var node = nodes[i];
                ib = context.createIndexBuffer(new Uint16Array(node.indices), usage, IndexDatatype.UNSIGNED_SHORT);
                va = context.createVertexArray(attributes, ib);
                vas.push({
                    id : node.coreId,
                    type : PrimitiveType.TRIANGLES,
                    vertexArray : va
                });
            }

            return vas;
        }
    }

    function getGeometryList(context, nodes) {
        var meshes = {};
        var length = nodes.length;
        for (var i = 0; i < length; ++i) {
            var node = nodes[i];
            if (node && node.type && node.type === 'geometry') {
                var geom = getGeometry(context, node);
                if (typeof geom !== 'undefined') {
                    for (var j = 0; j < geom.length; ++j) {
                        meshes[geom[j].id] = {
                                type : geom[j].type,
                                vertexArray : geom[j].vertexArray
                        };
                    }
                }
            }
        }
        return meshes;
    }

    Model.prototype.update = function(context, sceneState) {
        if (typeof this._sp === 'undefined') {
            var vs = '';
            vs += 'attribute vec3 position;';
            vs += 'attribute vec3 normal;';
            vs += 'uniform float u_scale;';
            vs += 'varying vec3 v_positionEC;';
            vs += 'varying vec3 v_normalEC;';
            vs += 'void main()';
            vs += '{';
            vs += '    v_positionEC = (agi_modelView * vec4(position, 1.0)).xyz;';
            vs += '    v_normalEC = normalize(agi_normal * normal);';
            vs += '    gl_Position = agi_modelViewProjection * vec4(u_scale * position, 1.0);';
            vs += '}';
            var fs = '';
            fs += 'varying vec3 v_positionEC;';
            fs += 'varying vec3 v_normalEC;';
            fs += 'void main()';
            fs += '{';
            fs += '    agi_material material;';
            fs += '    material.normal = normalize(v_normalEC);';
            fs += '    material.specular = 0.2;';
            fs += '    material.diffuse = vec3(0.34902, 0.0156863, 0.0823529);';
            fs += '    material.emission = vec3(0.0);';
            fs += '    material.alpha = 1.0;';
            fs += '    vec3 positionToEyeEC = normalize(-v_positionEC);';
            fs += '    gl_FragColor = agi_lightValuePhong(agi_sunDirectionEC, positionToEyeEC, material);';
            fs += '}';

            var attributeIndices = {
                    position : 0,
                    normal : 1
            };
            this._sp = context.getShaderCache().getShaderProgram(vs, fs, attributeIndices);

            this._rs = context.createRenderState({
                depthTest : {
                    enabled : true
                },
                cullFace : {
                    enabled : true
                }
            });
        }

        if (typeof this._vas === 'undefined') {
            var nodes = this._data.nodes;
            if (typeof nodes !== 'undefined') {
                var length = nodes.length;
                for (var i = 0; i < length; ++i) {
                    var node = nodes[i];
                    if (node && node.type && node.type === 'library') {
                        this._vas = getGeometryList(context, node.nodes);
                        break;
                    }
                }
            }
        }
    };

    Model.prototype.render = function(context) {
        if (typeof this._vas !== 'undefined') {
            for (var id in this._vas) {
                if (this._vas.hasOwnProperty(id)) {
                    var primitive = this._vas[id];
                    context.draw({
                        primitiveType : primitive.type,
                        shaderProgram : this._sp,
                        uniformMap    : this._drawUniforms,
                        vertexArray   : primitive.vertexArray,
                        renderState   : this._rs
                    });
                }
            }
        }
    };

    return Model;
});