import {TestType} from "./types";

function f(): TestType {
    return {
        num: 1,
        str: '',
    }
}

console.log(f());
